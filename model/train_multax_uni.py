import os
import torch
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, random_split
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import CSVLogger

import numpy as np

from fairseq.modules import LayerNorm
from fairseq.modules.quant_noise import quant_noise as apply_quant_noise_

import sys
from gly_encoder import GlyEncoder,GlyEncoder_el
from galy_dataset import GlyDataset_ml
from sklearn.metrics import f1_score

num_gpu=1
# num_gpu=[1]

encode_mode=1
batchsize=512
# batchsize=1
emb_dim=384
ffn_dim=1536
num_head=12
num_layer=12
vocab_dim=417
anomeric_dim=3
conn_dim=9
dense_dim=128
dense_dim_k=256
ddk0=4
ddk1=6
ddk2=28
ddk3=84
ddk4=206
ddk5=408
mask_prob=0.25
log_name='encode_1_s28_12_384_12_b512_ml_mlmamc_25_new_nb'
bimmu=False
save_path='/glycan/glycp_'+log_name
sycn=False
class glycanbert(pl.LightningModule):

    def __init__(self):
        super().__init__()
        if encode_mode==1:
            self.encoder = GlyEncoder(
                padding_idx=0,
                vocab_size=vocab_dim,
                num_branch = 2**11,
                max_branch=50,
                connections=2**9,
                aromatic=4,
                num_encoder_layers=num_layer,
                embedding_dim=emb_dim,
                ffn_embedding_dim=ffn_dim,
                num_attention_heads=num_head,
                dropout=0.1,
                attention_dropout=0.1,
                activation_dropout=0.0,
                layerdrop=0.1,
                max_seq_len=64,
                encoder_normalize_before=True,
                apply_bert_init=True,
                activation_fn="gelu"
            )
        elif encode_mode==2:
            self.encoder = GlyEncoder_el(
                padding_idx=0,
                vocab_size=vocab_dim,
                num_branch = 11,
                max_dep=50,
                connections=9,
                aromatic=3,
                parent=10,
                num_encoder_layers=num_layer,
                embedding_dim=emb_dim,
                ffn_embedding_dim=ffn_dim,
                num_attention_heads=num_head,
                dropout=0.1,
                attention_dropout=0.1,
                activation_dropout=0.0,
                layerdrop=0.1,
                max_seq_len=64,
                encoder_normalize_before=True,
                apply_bert_init=True,
                activation_fn="gelu"
            )
            
    
        self.activation_fn = nn.GELU()
        
#         ml head
        self.mlmproj0=torch.nn.Conv1d(emb_dim, emb_dim, 1)
        self.mlmproj1=torch.nn.Conv1d(emb_dim, vocab_dim, 1)
        self.layer_norm = LayerNorm(emb_dim)
        
        #         ma head
        self.mamproj0=torch.nn.Conv1d(emb_dim, emb_dim, 1)
        self.mamproj1=torch.nn.Conv1d(emb_dim, anomeric_dim, 1)
        
        #         mc head
        self.mcmproj0=torch.nn.Conv1d(emb_dim, emb_dim, 1)
        self.mcmproj1=torch.nn.Conv1d(emb_dim, conn_dim, 1)
        
#         cls head
        self.denselink = nn.Linear(emb_dim, dense_dim)
        self.densekingdom = nn.Linear(emb_dim, dense_dim_k)
        self.denseimmu = nn.Linear(emb_dim, dense_dim)
        
        self.dropout = nn.Dropout(p=0.0)
        self.out_link = nn.Linear(dense_dim, 3)
        if bimmu:
            self.out_immu = nn.Linear(dense_dim, 1)
        else:
            self.out_immu = nn.Linear(dense_dim, 2)
        
        self.out_k0 = nn.Linear(dense_dim_k, ddk0)
        self.out_k1 = nn.Linear(dense_dim_k, ddk1)
        self.out_k2 = nn.Linear(dense_dim_k, ddk2)
        self.out_k3 = nn.Linear(dense_dim_k, ddk3)
        self.out_k4 = nn.Linear(dense_dim_k, ddk4)
        self.out_k5 = nn.Linear(dense_dim_k, ddk5)


    def forward(self, x):
        # in lightning, forward defines the prediction/inference actions
        embedding,rep = self.encoder(x)
        
#       proj back
        emb=embedding[0].transpose(0,1)
        emb=emb.transpose(1,2)
        pb = self.mlmproj0(emb)
        pb = self.activation_fn(pb)
        pb=pb.transpose(1,2)
        pb = self.layer_norm(pb)
        # project back to size of vocabulary with bias
        pb=pb.transpose(1,2)
        pb = self.mlmproj1(pb) 
        
        pba = self.mamproj0(emb)
        pba = self.activation_fn(pba)
        pba=pba.transpose(1,2)
        pba = self.layer_norm(pba)
        # project back to size of vocabulary with bias
        pba=pba.transpose(1,2)
        pba = self.mamproj1(pba) 
        
        pbc = self.mcmproj0(emb)
        pbc = self.activation_fn(pbc)
        pbc=pbc.transpose(1,2)
        pbc = self.layer_norm(pbc)
        # project back to size of vocabulary with bias
        pbc=pbc.transpose(1,2)
        pbc = self.mcmproj1(pbc) 
        
#         cls
        
        cl = self.dropout(rep)
        cl = self.denselink(cl)
        cl = self.activation_fn(cl)
        cl = self.dropout(cl)
        cl = self.out_link(cl)
        
        ck = self.dropout(rep)
        ck = self.densekingdom(ck)
        ck = self.activation_fn(ck)
        ck = self.dropout(ck)
        ck0 = self.out_k0(ck)
        ck1 = self.out_k1(ck)
        ck2 = self.out_k2(ck)
        ck3 = self.out_k3(ck)
        ck4 = self.out_k4(ck)
        ck5 = self.out_k5(ck)
        
        ci = self.dropout(rep)
        ci = self.denseimmu(ci)
        ci = self.activation_fn(ci)
        ci = self.dropout(ci)
        ci = self.out_immu(ci)
        return pb,pba,pbc,cl,ci,ck0,ck1,ck2,ck3,ck4,ck5
    

    def training_step(self, batch, batch_idx):
        # training_step defined the train loop. It is independent of forward
        
        x, y,tl0,tl1,tl2,tl3,tl4,tl5,mcbl = batch
        
        tl0=tl0.squeeze()
        tl1=tl1.squeeze()
        tl2=tl2.squeeze()
        tl3=tl3.squeeze()
        tl4=tl4.squeeze()
        tl5=tl5.squeeze()
        
        bs=x.size()[0]
        
        mlmy=x[:,:,0].clone()
        
        mask=torch.full(mlmy.size(),False)
        for i in range(0,bs):
            ordlen=torch.sum(mlmy[i,:]>3).cpu()
            num_mask=max(int(ordlen*mask_prob+np.random.rand()),1)
            mask_idc = np.random.choice(np.arange(1,ordlen), num_mask, replace=False)
            for m in mask_idc:
                mask[i,m]=True
                p=np.random.rand()
                if p<=0.1:
                    x[i,m,0]=np.random.randint(4,417)
                elif p<=0.9:
                    x[i,m,0]=1

        mamy=x[:,:,1].clone()
        maska=torch.full(mamy.size(),False)
        for i in range(0,bs):
            ordalen=torch.sum(mamy[i,:]>0).cpu()
            ordlen=torch.sum(mlmy[i,:]>3).cpu()
            num_mask=int(ordalen*mask_prob+np.random.rand())
            if len([j for j in np.arange(1,ordlen) if x[i,j,1]>0])<num_mask:
                continue
            mask_idc = np.random.choice([j for j in np.arange(1,ordlen) if x[i,j,1]>0], num_mask, replace=False)
            for m in mask_idc:
                maska[i,m]=True
                p=np.random.rand()
                if p<=0.1:
                    x[i,m,1]=np.random.randint(0,3)
                elif p<=0.9:
                    x[i,m,1]=3

        mcmy=x[:,:,2].clone()
        maskc=torch.full(mcmy.size(),False)
        for i in range(0,bs):
            ordclen=torch.sum(mcmy[i,:]>0).cpu()
            ordlen=torch.sum(mlmy[i,:]>3).cpu()
            num_mask=int(ordclen*mask_prob+np.random.rand())
            if len([j for j in np.arange(1,ordlen) if x[i,j,2]>0])<num_mask:
                continue
            mask_idc = np.random.choice([j for j in np.arange(1,ordlen) if x[i,j,2]>0], num_mask, replace=False)
            for m in mask_idc:
                maska[i,m]=True
                p=np.random.rand()
                if p<=0.1:
                    x[i,m,2]=np.random.randint(1,500)
                elif p<=0.9:
                    x[i,m,2]=511
        
        lml,lal,lmc,lcl,lci,lck0,lck1,lck2,lck3,lck4,lck5=self(x)
        
        lml=lml.transpose(1,2)
        loss = F.cross_entropy(lml[mask],mlmy[mask])
        self.log('mlm', loss, sync_dist=sycn)

        loss2=0
        lal=lal.transpose(1,2)
        loss2 += F.cross_entropy(lal[maska],mamy[maska])

        lmc=lmc.transpose(1,2)
        loss2+=F.binary_cross_entropy_with_logits(lmc[maskc,:],mcbl[maskc,:].type_as(lmc),pos_weight=(torch.numel(mcbl[maskc,:])-mcbl[maskc,:].sum())/mcbl[maskc,:].sum())
        
        
        lw=torch.tensor([1/59,1/721,1/509]).type_as(lcl)
        yl=y[:,0].clone()
        syl=yl>-1
        if syl.any():
            loss2+= F.cross_entropy(lcl[syl],yl[syl],weight=lw)
            accuracy = (torch.softmax(lcl[syl], dim=1).argmax(dim=1) == yl[syl]).sum().float() / float( yl[syl].size(0) )
            self.log('link_a', accuracy, sync_dist=sycn)

        yi=y[:,1].clone()
        syi=yi>-1
        if syi.any():
            if bimmu:
                accuracy = ((torch.sigmoid(lci[syi,:])>0.5).view(-1)==yi[syi].view(-1)).sum().float()/float(yi[syi].size(0))
                loss2+=F.binary_cross_entropy_with_logits(lci[syi,:].view(-1),yi[syi].type_as(lci),pos_weight=(torch.numel(yi[syi])-yi[syi].sum())/yi[syi].sum())
            else:
                accuracy = (torch.softmax(lci[syi], dim=1).argmax(dim=1) == yi[syi]).sum().float() / float( yi[syi].size(0) )
                loss2+= F.cross_entropy(lci[syi],yi[syi])
            self.log('immu_a', accuracy, sync_dist=sycn)

        syk0=(torch.sum(tl0,dim=1)>0)
        if syk0.any():
            loss2+= F.binary_cross_entropy_with_logits(lck0[syk0,:].squeeze(),tl0[syk0,:].squeeze(),pos_weight=(torch.numel(tl0[syk0,:])-tl0[syk0,:].sum())/tl0[syk0,:].sum())
            accuracy =f1_score(tl0[syk0].cpu().numpy(),torch.gt(torch.sigmoid(lck0[syk0]),0.5).detach().cpu().numpy(),average='samples')
            self.log('t0_f', accuracy, sync_dist=sycn)
            
        syk1=(torch.sum(tl1,dim=1)>0)
        if syk1.any():
            loss2+= F.binary_cross_entropy_with_logits(lck1[syk1,:].squeeze(),tl1[syk1,:].squeeze(),pos_weight=(torch.numel(tl1[syk1,:])-tl1[syk1,:].sum())/tl1[syk1,:].sum())
            accuracy =f1_score(tl1[syk1].cpu().numpy(),torch.gt(torch.sigmoid(lck1[syk1]),0.5).detach().cpu().numpy(),average='samples')
            self.log('t1_f', accuracy, sync_dist=sycn)
            
        syk2=(torch.sum(tl2,dim=1)>0)
        if syk2.any():
            loss2+= F.binary_cross_entropy_with_logits(lck2[syk2,:].squeeze(),tl2[syk2,:].squeeze(),pos_weight=(torch.numel(tl2[syk2,:])-tl2[syk2,:].sum())/tl2[syk2,:].sum())
            accuracy =f1_score(tl2[syk2].cpu().numpy(),torch.gt(torch.sigmoid(lck2[syk2]),0.5).detach().cpu().numpy(),average='samples')
            self.log('t2_f', accuracy, sync_dist=sycn)
            
        syk3=(torch.sum(tl3,dim=1)>0)
        if syk3.any():
            loss2+= F.binary_cross_entropy_with_logits(lck3[syk3,:].squeeze(),tl3[syk3,:].squeeze(),pos_weight=(torch.numel(tl3[syk3,:])-tl3[syk3,:].sum())/tl3[syk3,:].sum())
            
        syk4=(torch.sum(tl4,dim=1)>0)
        if syk4.any():
            loss2+= F.binary_cross_entropy_with_logits(lck4[syk4,:].squeeze(),tl4[syk4,:].squeeze(),pos_weight=(torch.numel(tl4[syk4,:])-tl4[syk4,:].sum())/tl4[syk4,:].sum())

        syk5=(torch.sum(tl5,dim=1)>0)
        if syk5.any():
            loss2+= F.binary_cross_entropy_with_logits(lck5[syk5,:].squeeze(),tl5[syk5,:].squeeze(),pos_weight=(torch.numel(tl5[syk5,:])-tl5[syk5,:].sum())/tl5[syk5,:].sum())
        
            
            
            
        self.log('total', loss+loss2, sync_dist=sycn)

        return loss+loss2
    
    def validation_step(self, batch, batch_idx):
        
        x, y,tl0,tl1,tl2,tl3,tl4,tl5,mcbl = batch
        
        tl0=tl0.squeeze()
        tl1=tl1.squeeze()
        tl2=tl2.squeeze()
        tl3=tl3.squeeze()
        tl4=tl4.squeeze()
        tl5=tl5.squeeze()
        
        bs=x.size()[0]
        
        mlmy=x[:,:,0].clone()
        
        mask=torch.full(mlmy.size(),False)
        for i in range(0,bs):
            ordlen=torch.sum(mlmy[i,:]>3).cpu()
            num_mask=max(int(ordlen*mask_prob+np.random.rand()),1)
            mask_idc = np.random.choice(np.arange(1,ordlen), num_mask, replace=False)
            for m in mask_idc:
                mask[i,m]=True
                p=np.random.rand()
                if p<=0.1:
                    x[i,m,0]=np.random.randint(0,417)
                elif p<=0.9:
                    x[i,m,0]=1
        
        lml,lal,lmc,lcl,lci,lck0,lck1,lck2,lck3,lck4,lck5=self(x)
        
        lml=lml.transpose(1,2)
        loss = F.cross_entropy(lml[mask],mlmy[mask])
        self.log('val_mlm', loss)
        accuracy = (torch.softmax(lml[mask], dim=1).argmax(dim=1) == mlmy[mask]).sum().float() / float( mlmy[mask].size(0) )
        self.log('val_mlm_a', accuracy, sync_dist=sycn)
        
        loss2=0
        lw=torch.tensor([1/59,1/721,1/509]).type_as(lcl)
        yl=y[:,0].clone()
        syl=yl>-1
        if syl.any():
            loss2+= F.cross_entropy(lcl[syl],yl[syl],weight=lw)
            accuracy = (torch.softmax(lcl[syl], dim=1).argmax(dim=1) == yl[syl]).sum().float() / float( yl[syl].size(0) )
            self.log('val_link_a', accuracy, sync_dist=sycn)
            
            
        yi=y[:,1].clone()
        syi=yi>-1
        if syi.any():
            if bimmu:
                accuracy = ((torch.sigmoid(lci[syi,:])>0.5).view(-1)==yi[syi].view(-1)).sum().float()/float(yi[syi].size(0))
                loss2+=F.binary_cross_entropy_with_logits(lci[syi,:].view(-1),yi[syi].type_as(lci),pos_weight=(torch.numel(yi[syi])-yi[syi].sum())/yi[syi].sum())
            else:
                accuracy = (torch.softmax(lci[syi], dim=1).argmax(dim=1) == yi[syi]).sum().float() / float( yi[syi].size(0) )
                loss2+= F.cross_entropy(lci[syi],yi[syi])
                
            self.log('val_immu_a', accuracy, sync_dist=sycn)

        syk0=(torch.sum(tl0,dim=1)>0)
        if syk0.any():
            loss2+= F.binary_cross_entropy_with_logits(lck0[syk0,:].squeeze(),tl0[syk0,:].squeeze(),pos_weight=(torch.numel(tl0[syk0,:])-tl0[syk0,:].sum())/tl0[syk0,:].sum())
            accuracy =f1_score(tl0[syk0].cpu().numpy(),torch.gt(torch.sigmoid(lck0[syk0]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t0_f', accuracy, sync_dist=sycn)
            
        syk1=(torch.sum(tl1,dim=1)>0)
        if syk1.any():
            loss2+= F.binary_cross_entropy_with_logits(lck1[syk1,:].squeeze(),tl1[syk1,:].squeeze(),pos_weight=(torch.numel(tl1[syk1,:])-tl1[syk1,:].sum())/tl1[syk1,:].sum())
            accuracy =f1_score(tl1[syk1].cpu().numpy(),torch.gt(torch.sigmoid(lck1[syk1]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t1_f', accuracy, sync_dist=sycn)
            
        syk2=(torch.sum(tl2,dim=1)>0)
        if syk2.any():
            loss2+= F.binary_cross_entropy_with_logits(lck2[syk2,:].squeeze(),tl2[syk2,:].squeeze(),pos_weight=(torch.numel(tl2[syk2,:])-tl2[syk2,:].sum())/tl2[syk2,:].sum())
            accuracy =f1_score(tl2[syk2].cpu().numpy(),torch.gt(torch.sigmoid(lck2[syk2]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t2_f', accuracy, sync_dist=sycn)
            
        syk3=(torch.sum(tl3,dim=1)>0)
        if syk3.any():
            loss2+= F.binary_cross_entropy_with_logits(lck3[syk3,:].squeeze(),tl3[syk3,:].squeeze(),pos_weight=(torch.numel(tl3[syk3,:])-tl3[syk3,:].sum())/tl3[syk3,:].sum())
            accuracy =f1_score(tl3[syk3].cpu().numpy(),torch.gt(torch.sigmoid(lck3[syk3]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t3_f', accuracy, sync_dist=sycn)
            
        syk4=(torch.sum(tl4,dim=1)>0)
        if syk4.any():
            loss2+= F.binary_cross_entropy_with_logits(lck4[syk4,:].squeeze(),tl4[syk4,:].squeeze(),pos_weight=(torch.numel(tl4[syk4,:])-tl4[syk4,:].sum())/tl4[syk4,:].sum())
            accuracy =f1_score(tl4[syk4].cpu().numpy(),torch.gt(torch.sigmoid(lck4[syk4]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t4_f', accuracy, sync_dist=sycn)
            
        syk5=(torch.sum(tl5,dim=1)>0)
        if syk5.any():
            loss2+= F.binary_cross_entropy_with_logits(lck5[syk5,:].squeeze(),tl5[syk5,:].squeeze(),pos_weight=(torch.numel(tl5[syk5,:])-tl5[syk5,:].sum())/tl5[syk5,:].sum())
            accuracy =f1_score(tl5[syk5].cpu().numpy(),torch.gt(torch.sigmoid(lck5[syk5]),0.5).detach().cpu().numpy(),average='samples')
            self.log('val_t5_f', accuracy, sync_dist=sycn)
            
        self.log('val_total', loss2+loss)
        

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(
            self.parameters(),
            betas=(0.9, 0.98),
            lr=5e-4,
            eps=1e-06, 
            weight_decay=0.01
        )
        return optimizer
    
    def train_dataloader(self):
        if encode_mode==1:
            train=GlyDataset_ml('train_28s_ml_uni.pkl')
        elif encode_mode==2:
            train=GlyDataset('train_28s_ne0.pkl')
        return DataLoader(train,batch_size=batchsize,collate_fn=train.collater,num_workers=4)
        
    def val_dataloader(self):
        if encode_mode==1:
            valid=GlyDataset_ml('valid_28s_ml_uni.pkl')
        elif encode_mode==2:
            valid=GlyDataset('valid_28s_ne0.pkl')
        return DataLoader(valid,batch_size=batchsize,collate_fn=valid.collater,num_workers=4)

gbt = glycanbert()

checkpoint_val = ModelCheckpoint(
    monitor='val_total',
    dirpath=save_path,
    filename='glybert-val-{epoch:02d}-{val_total:.3f}',
#     period=1,
    save_top_k=3,
    mode='min',
)

checkpoint_train = ModelCheckpoint(
    monitor='total',
    dirpath=save_path,
    filename='glybert-train-{epoch:02d}-{total:.3f}',
#     period=1,
    save_top_k=3,
    mode='min',
)

logger = CSVLogger("logs", name=log_name)

trainer = pl.Trainer(
    max_epochs=200,
    gpus=num_gpu, 
    precision=16,
    accumulate_grad_batches=1,
    gradient_clip_val=0,
    check_val_every_n_epoch=5,
    callbacks=[checkpoint_train,checkpoint_val],
    logger=logger
)
trainer.fit(gbt)