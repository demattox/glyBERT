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
sys.path.append('/home/bowen/glycan/galy')
from gly_encoder import GlyEncoder,GlyEncoder_el
from galy_dataset import GlyDataset_cfg_cv
from sklearn.metrics import f1_score,roc_curve,auc,roc_auc_score

num_gpu=1
# num_gpu=[1]

encode_mode=1
batchsize=128
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
mask_prob=0.35

densetune0=256
densetune1=128
densetuneout=1
file=sys.argv[1]
logfolder='cfg3_logs'
log_name=file
save_path='/home/bowen/glycan/cfgcv3/'+log_name
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
        self.out_immu = nn.Linear(dense_dim, 2)
        
        self.out_k0 = nn.Linear(dense_dim_k, ddk0)
        self.out_k1 = nn.Linear(dense_dim_k, ddk1)
        self.out_k2 = nn.Linear(dense_dim_k, ddk2)
        self.out_k3 = nn.Linear(dense_dim_k, ddk3)
        self.out_k4 = nn.Linear(dense_dim_k, ddk4)
        self.out_k5 = nn.Linear(dense_dim_k, ddk5)
        
        cp=torch.load('checkpoint/glycp_encode_1_s28_12_384_12_b512_ml_mlmamc_15_new_nb/glybert-val-epoch=184-val_total=2.986.ckpt')
        self.load_state_dict(cp['state_dict'])
        
        self.tune0 = nn.Linear(emb_dim, densetune0)
        self.tune1 = nn.Linear(densetune0, densetune1)
        self.tune2 = nn.Linear(densetune1, densetuneout)

    def forward(self, x):
        # in lightning, forward defines the prediction/inference actions
        embedding,rep = self.encoder(x)
    
        cl = self.dropout(rep)
        cl = self.tune0(cl)
        cl = self.activation_fn(cl)
        cl = self.dropout(cl)
        cl = self.tune1(cl)
        cl = self.activation_fn(cl)
        cl = self.dropout(cl)
        cl = self.tune2(cl)
        
        return cl
        
        

    

    def training_step(self, batch, batch_idx):

        
        x, y = batch
        
        
        
        cy=self(x)
        
        
        loss=F.binary_cross_entropy_with_logits(cy,y.type_as(cy).unsqueeze(1),pos_weight=(torch.numel(y)-y.sum())/y.sum())
        try:
            accuracy =roc_auc_score(y.cpu().numpy(),torch.sigmoid(cy).detach().cpu().numpy())
            self.log('auc', accuracy, sync_dist=sycn)
        except:
            pass
            
        self.log('total', loss, sync_dist=sycn)

        return loss
    
    def validation_step(self, batch, batch_idx):
        
        x, y = batch
        
        cy=self(x)
        
        
        loss=F.binary_cross_entropy_with_logits(cy,y.type_as(cy).unsqueeze(1),pos_weight=(torch.numel(y)-y.sum())/y.sum())

        try:
            accuracy =roc_auc_score(y.cpu().numpy(),torch.sigmoid(cy).detach().cpu().numpy())
            self.log('val_auc', accuracy, sync_dist=sycn)
        except:
            pass
            
        self.log('val_total', loss)
        

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
            train=GlyDataset_cfg_cv('cfg_cv2/train_'+file+'.pkl')

        return DataLoader(train,batch_size=batchsize,collate_fn=train.collater,num_workers=2)
        
    def val_dataloader(self):
        if encode_mode==1:
            valid=GlyDataset_cfg_cv('cfg_cv2/test_'+file+'.pkl')
        
        return DataLoader(valid,batch_size=batchsize,collate_fn=valid.collater,num_workers=2)
    


gbt = glycanbert()

checkpoint_val = ModelCheckpoint(
    monitor='val_total',
    dirpath=save_path,
    filename='glybert-val-{epoch:02d}-{val_total:.3f}',
    save_top_k=3,
    mode='min',
)

checkpoint_train = ModelCheckpoint(
    monitor='total',
    dirpath=save_path,
    filename='glybert-train-{epoch:02d}-{total:.3f}',
    save_top_k=3,
    mode='min',
)

logger = CSVLogger(logfolder, name=log_name)


trainer = pl.Trainer(
    max_epochs=100,
    gpus=num_gpu, 
    precision=16,
    accumulate_grad_batches=1,
    gradient_clip_val=0,
    check_val_every_n_epoch=5,
    logger=logger
)

trainer.fit(gbt)