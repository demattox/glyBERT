from glycanbert_mlmamc import glycanbert
from torch.utils.data import DataLoader
import numpy as np
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
import collections
from sklearn.manifold import TSNE
from galy_dataset import GlyDataset
# from gly_encoder import GlyEncoder,GlyEncoder_el
from galy_dataset import GlyDataset_ml
from torch.distributions.categorical import Categorical
model=glycanbert.load_from_checkpoint('/checkpoint/glycp_encode_1_s28_12_384_12_b512_ml_mlmamc_15_new_nb/glybert-val-epoch=184-val_total=2.986.ckpt')
model.eval()

model3=glycanbert.load_from_checkpoint('/checkpoint/glycp_encode_1_s28_12_384_12_b512_ml_mlmamc_15_new_nb/glybert-val-epoch=184-val_total=2.986.ckpt')
model3.eval()

import pickle
import torch.nn.functional as F

class Identity(torch.nn.Module):
    def __init__(self):
        super(Identity, self).__init__()
        
    def forward(self, x):
        return x
    
class Zero(torch.nn.Module):
    def __init__(self):
        super(Zero, self).__init__()
        
    def forward(self, x):
        return 0

def enc2x(start):
    allx=start.long()
    padding_mask = allx[:,:,0].eq(model.encoder.padding_idx)
    if not model.encoder.traceable and not padding_mask.any():
        padding_mask = None

#     print(padding_mask)
    x = model.encoder.embed_tokens(allx[:,:,0])

    if model.encoder.embed_scale is not None:
        x = x * model.encoder.embed_scale

    if model.encoder.connection_embeddings is not None:
        x = x + model.encoder.connection_embeddings(allx[:,:,2])

    if model.encoder.branch_embeddings is not None:
        x = x + model.encoder.branch_embeddings(allx[:,:,3])

    if model.encoder.depth_embeddings is not None:
        x = x + model.encoder.depth_embeddings(allx[:,:,5])

    if model.encoder.aromatic_embeddings is not None:
        x = x + model.encoder.aromatic_embeddings(allx[:,:,1])

    if model.encoder.parent_embeddings is not None:
        x = x + model.encoder.parent_embeddings(allx[:,:,4])

    if model.encoder.quant_noise is not None:
        x = model.encoder.quant_noise(x)

    return x

def getgrad(start,t,ml=0,target=None,p=0):
    x=enc2x(start)
    nx=x.data
    nx.requires_grad=True
    lml,lal,lmc,lcl,lci,lck0,lck1,lck2,lck3,lck4,lck5,emb,rep=model3(nx)
    model3.zero_grad()
    if ml==0:
        loss=F.cross_entropy(lci,torch.FloatTensor([t]).long())
    elif ml==1:
        loss=F.mse_loss(rep,target)
    loss.backward(retain_graph=True)
    d=nx.grad.data
    if p==1:
        print(lci)
    return lml,lal,nx,d

def update(nx,d,a,t,ml=0,target=None,p=0):
    nx=nx.data-d*a
    nx.requires_grad=True
    lml,lal,lmc,lcl,lci,lck0,lck1,lck2,lck3,lck4,lck5,emb,rep=model3(nx)
    model3.zero_grad()
    if ml==0:
        loss=F.cross_entropy(lci,torch.FloatTensor([t]).long())
    elif ml==1:
        loss=F.mse_loss(rep,target)
    loss.backward(retain_graph=True)
    d=nx.grad.data
    if p==1:
        print(lci)
    return lml,lal,nx,d

    
model3.encoder.gen=True

model3.encoder.embed_tokens=Identity()


model3.encoder.connection_embeddings=Zero()
model3.encoder.branch_embeddings=Zero()
model3.encoder.depth_embeddings=Zero()
model3.encoder.aromatic_embeddings =Zero()
model3.encoder.parent_embeddings=Zero()
rnum=100
rl=10
shortcutoff=9
with open('alignstruct.pkl', 'rb') as inFH:
    datanewall = pickle.load(inFH)
    
qualigly=[2,3,4,8,10,13,22]
for sele in qualigly:
    nc=[rl]*len(datanew[sele]['No'])
    trec=[]
    pathrec=[]
    pathall=[]
    see=0
    with open('../data/processed/glycan_common_16048_may22.pickle', 'rb') as inFH:
        dataall = pickle.load(inFH)
    
    currdict={}
    currdict['Yes']=[]
    currdict['No']=[]

    for key in ['No','Yes']:
        for ids in datanew[sele][key]:
            rep=dataall[ids][0]
            rep[:,0]+=2
            mask=rep[:,0]>3

            currdict[key].append(rep[mask,0].tolist())
    
    with open('../data/processed/glycan_common_16048_may22.pickle', 'rb') as inFH:
        dataall = pickle.load(inFH)
    alldict=[]
    for ids in datanewall[sele]:
        rep=dataall[ids][0]
        rep[:,0]+=2
        mask=rep[:,0]>3

        alldict.append(rep[mask,0].tolist())
    
    assert currdict['Yes'][0] in alldict
                
    with open('../data/processed/glycan_common_16048_may22.pickle', 'rb') as inFH:
        dataall = pickle.load(inFH)

    for j,key in enumerate(datanew[sele]['No']):     
        # start no
        start=torch.from_numpy(dataall[key][0]).unsqueeze(0)
        start[:,:,0]+=2
        t=1
        lr=0.1
        mask=start[:,:,0]<=3
        print(start[:,:,0][start[:,:,0]>3])
        it=rl
        shortpath=[]
        restpath=[]
        for k in range(rnum):
            temppath=[]
            temppath.append(start[:,:,0][start[:,:,0]>3].tolist())
            for i in range(rl):
                if i==0:
                    cstart=start.clone()
                else:
                    cstart[0,pick,1]=3
                    cstart[0,pick,2]=511
                lml0,lal0,nx,d=getgrad(cstart,t)
                p0=F.softmax(lml0,dim=1)
                p0m=p0[0,cstart[0,:,0],range(p0.size(2))]

                lml1,lal1,nx,d=update(nx,d,lr,t)

                p1=F.softmax(lml1,dim=1)

                p1m=p1[0,cstart[0,:,0],range(p1.size(2))]

                pdiff=p0m-p1m
                pdiff[mask.squeeze()]=0
                pdiff=torch.abs(pdiff)
                normpdiff=pdiff/pdiff.sum()
                p=Categorical(normpdiff.squeeze())

                pick=p.sample()
                pd=p1[0,:,pick]
                pd[pick]=0
                m=Categorical(pd)

                new=m.sample()

                while new==cstart[0,pick,0] or new <=3:
                    new=m.sample()

                cstart[0,pick,0]=new
                temppath.append(cstart[:,:,0][~mask].tolist())
            
                if cstart[:,:,0][~mask].tolist() in alldict:
                    see+=1
                if cstart[:,:,0][~mask].tolist() in currdict['Yes']:
                    it=min(it,i+1)
                    break
            trec.append(i+1)
            
            if len(temppath)<shortcutoff:
                if temppath not in shortpath:
                    shortpath.append(temppath)
            else:
                if temppath not in restpath:
                    restpath.append(temppath)
                    
        nc[j]=min(nc[j],it)
        pathrec.append(shortpath)
        pathall.append(restpath)
        
    pickle.dump([nc,trec,pathrec,pathall,see],open('genrecords9/'+str(sele)+'.pkl','wb'))
