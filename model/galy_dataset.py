import os
import json
import numpy as np

import torch
import torch.nn.functional as F
import pickle

from fairseq.data import FairseqDataset, data_utils
from torch.utils.data.dataloader import default_collate
import sys


class GlyDataset(FairseqDataset):
    def __init__(self, file,shuffle=False):
        self.file = file
        with open(self.file, 'rb') as inFH:
            self.data = pickle.load(inFH)
        self.keys=list(self.data.keys())
        
    def __getitem__(self, index):
        
        return torch.from_numpy(self.data[self.keys[index]][0]).long(),torch.from_numpy(self.data[self.keys[index]][1])

    def __len__(self):
        return len(self.data)

    def num_tokens(self, index):
        return self.data[self.keys[index]][0].shape[0]
    
    def size(self, index):
        return data[self.keys[index]][0].shape

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
#         num_objects = [features.shape[0] for features, _ in samples]
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        label_samples_padded = []

        for (features,label), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            label_samples_padded.append(label)

        return default_collate(feature_samples_padded),default_collate(label_samples_padded)

    
    
class GlyDataset_ml(FairseqDataset):
    def __init__(self, file,shuffle=False):
        self.file = file
        with open(self.file, 'rb') as inFH:
            self.data = pickle.load(inFH)
        self.keys=list(self.data.keys())

    def __getitem__(self, index):
        tdata=torch.from_numpy(self.data[self.keys[index]][0]).long()
        lilabel=torch.from_numpy(self.data[self.keys[index]][1])
        lk0=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][0] if i>-1]).unsqueeze(0)
        lk1=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][1] if i>-1]).unsqueeze(0)
        lk2=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][2] if i>-1]).unsqueeze(0)
        lk3=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][3] if i>-1]).unsqueeze(0)
        lk4=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][4] if i>-1]).unsqueeze(0)
        lk5=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][2][5] if i>-1]).unsqueeze(0)
        
        tlk0 = torch.zeros(lk0.size(0), 4).scatter_(1, lk0.type(torch.int64), 1.)
        tlk1 = torch.zeros(lk1.size(0), 6).scatter_(1, lk1.type(torch.int64), 1.)
        tlk2 = torch.zeros(lk2.size(0), 28).scatter_(1, lk2.type(torch.int64), 1.)
        tlk3 = torch.zeros(lk3.size(0), 84).scatter_(1, lk3.type(torch.int64), 1.)
        tlk4 = torch.zeros(lk4.size(0), 206).scatter_(1, lk4.type(torch.int64), 1.)
        tlk5 = torch.zeros(lk5.size(0), 408).scatter_(1, lk5.type(torch.int64), 1.)
        
        binary_repr_v = np.vectorize(np.binary_repr)
        cl=torch.from_numpy(np.array([list(b) for b in binary_repr_v(tdata[:,2],9)]).astype(np.int8))
        
        return tdata,lilabel,tlk0,tlk1,tlk2,tlk3,tlk4,tlk5,cl

    def __len__(self):
        return len(self.data)

    def num_tokens(self, index):
        return self.data[self.keys[index]][0].shape[0]
    
    def size(self, index):
        return data[self.keys[index]][0].shape

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
        num_objects = [features.shape[0] for features,_,_,_,_,_,_,_,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        label_samples_padded = []
        l0=[]
        l1=[]
        l2=[]
        l3=[]
        l4=[]
        l5=[]
        cl=[]

        for (features,label,ll0,ll1,ll2,ll3,ll4,ll5,cll), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            label_samples_padded.append(label)
            l0.append(ll0)
            l1.append(ll1)
            l2.append(ll2)
            l3.append(ll3)
            l4.append(ll4)
            l5.append(ll5)
            pcll= F.pad(cll, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            cl.append(pcll)
            
        return default_collate(feature_samples_padded),default_collate(label_samples_padded),default_collate(l0),default_collate(l1),default_collate(l2),default_collate(l3),default_collate(l4),default_collate(l5),default_collate(cl)
    
    
class GlyDataset_cfg(FairseqDataset):
    def __init__(self, file,shuffle=False):
        self.file = file
        
        with open('encoding_pickles/CFG_glycans_535_feb9.pickle', 'rb') as inFH:
            self.data = pickle.load(inFH)
            
        with open(self.file, 'rb') as inFH:
            self.label = pickle.load(inFH)
        
        
        self.keys=list(self.label.keys())


    def __getitem__(self, index):
        
        return torch.from_numpy(self.data[self.keys[index]][0]).long(),torch.from_numpy(self.label[self.keys[index]])

    def __len__(self):
        return len(self.label)

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        label_samples_padded = []

        for (features,label), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            label_samples_padded.append(label)

        return default_collate(feature_samples_padded),default_collate(label_samples_padded)
    
    
class GlyDataset_cfg_cv(FairseqDataset):
    def __init__(self, file,shuffle=False):
        self.file = file
        
#         with open('encoding_pickles/CFG_glycans_535_feb9.pickle', 'rb') as inFH:
        with open('encoding_pickles/CFG_glycans_535_apr12.pickle', 'rb') as inFH:
            self.data = pickle.load(inFH)
            
        with open(self.file, 'rb') as inFH:
            self.label = pickle.load(inFH)
        
        
        self.keys=list(self.label.keys())


    def __getitem__(self, index):
        
        return torch.from_numpy(self.data[self.keys[index]][0]).long(),torch.from_numpy(np.asarray(self.label[self.keys[index]]))

    def __len__(self):
        return len(self.label)

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
#         num_objects = [features.shape[0] for features, _ in samples]
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        label_samples_padded = []

        for (features,label), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            label_samples_padded.append(label)

        return default_collate(feature_samples_padded),default_collate(label_samples_padded)
    
class GlyDataset_tax_tune(FairseqDataset):
    def __init__(self, file,level=0,shuffle=False):
        self.file = file
        with open(self.file, 'rb') as inFH:
            self.data = pickle.load(inFH)
        self.keys=list(self.data.keys())
        self.l=level

    def __getitem__(self, index):
        tdata=torch.from_numpy(self.data[self.keys[index]][0]).long()
        lk0=torch.FloatTensor([int(i) for i in self.data[self.keys[index]][1] if i>-1]).unsqueeze(0)
        
        ls=[4,3,28,84,206,408]
        
        tlk0 = torch.zeros(lk0.size(0), ls[self.l]).scatter_(1, lk0.type(torch.int64), 1.)
        
        return tdata,tlk0

    def __len__(self):
        return len(self.data)

    def num_tokens(self, index):
        return self.data[self.keys[index]][0].shape[0]
    
    def size(self, index):
        return data[self.keys[index]][0].shape

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        l0=[]
        

        for (features,ll0), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            l0.append(ll0)
            
            
        return default_collate(feature_samples_padded),default_collate(l0)
    
    
class GlyDataset_tax_tune_single(FairseqDataset):
    def __init__(self, file,level=0,shuffle=False):
        self.file = file
        with open(self.file, 'rb') as inFH:
            self.data = pickle.load(inFH)
        self.keys=list(self.data.keys())
        self.l=level

    def __getitem__(self, index):
        tdata=torch.from_numpy(self.data[self.keys[index]][0]).long()
        lk0=torch.FloatTensor([self.data[self.keys[index]][1][0]])
                        
        return tdata,lk0.long()

    def __len__(self):
        return len(self.data)

    def num_tokens(self, index):
        return self.data[self.keys[index]][0].shape[0]
    
    def size(self, index):
        return data[self.keys[index]][0].shape

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        l0=[]
        

        for (features,ll0), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            l0.append(ll0)
            
            
        return default_collate(feature_samples_padded),default_collate(l0)
    
    
class GlyDataset_tax_tune_single_jcgg(FairseqDataset):
    def __init__(self, file,level=0,shuffle=False):
        self.file = file
        with open(self.file, 'rb') as inFH:
            self.data = pickle.load(inFH)
        self.keys=list(self.data.keys())
        self.l=level

    def __getitem__(self, index):
        tdata=torch.from_numpy(self.data[self.keys[index]][0]).long()
        lk0=torch.FloatTensor([self.data[self.keys[index]][1]])
                        
        return tdata,lk0

    def __len__(self):
        return len(self.data)

    def num_tokens(self, index):
        return self.data[self.keys[index]][0].shape[0]
    
    def size(self, index):
        return data[self.keys[index]][0].shape

    def ordered_indices(self):
        if self.shuffle:
            indices = np.random.permutation(len(self))
        else:
            indices = np.arange(len(self))
            
        return indices

    def collater(self, samples):
        num_objects = [features.shape[0] for features,_ in samples]
        max_objects = max(num_objects)

        feature_samples_padded = []
        l0=[]
        

        for (features,ll0), n in zip(samples, num_objects):
            features_padded = F.pad(features, pad=[0,0,0, max_objects-n], mode='constant', value=0)
            feature_samples_padded.append(features_padded)
            l0.append(ll0)
            
            
        return default_collate(feature_samples_padded),default_collate(l0)