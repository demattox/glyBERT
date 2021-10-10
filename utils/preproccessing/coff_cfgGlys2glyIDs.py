#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:39:53 2021

@author: dmattox
"""

import re, os #, collections

homeDir = '/Users/dmattox/cbk/sugarTrees/'

coffDir = homeDir + 'glyPredict/coff_et_al_test-train_splits/'

cfgFile = homeDir + 'data/Clean_Lectin_intensities_67uni.txt'

embFile = homeDir + 'glyPredict/CFG_glycans_535_feb9_emb.csv'

embKeys = {}
with open(embFile, 'r') as inFH:
    for line in inFH:
        line = line.strip().split(',')
        if 'd0' not in line:
            k = line[0]
            embKeys[tuple([int(l) for l in k.split('_')])] = k

def getEmb(glyID, embDict = embKeys):
    return [v for k,v in embDict.items() if glyID in k][0]

getEmb(9)
getEmb(174)


cfgGlys = []
cfgGly2cfgID = dict()

with open(cfgFile, 'r') as inFH: # Read in sugar base csv file
    for line in inFH:
        line = line.strip()
        line = [l.strip() for l in line.split('\t')]
        if line[0] != 'glyID':
            gly = line[1].strip()
            gly = re.sub('[[:blank:]]', '', gly)
            gly = re.sub(' ', '', gly)
            gly = re.sub(r'\xa0', '', gly) # Non-breaking spaces for whatever dang reason
            cfgGlys.append(gly)
            cfgGly2cfgID[gly] = int(line[0])
            
coffGlys = []

trainSplitFiles = os.listdir(coffDir + 'single_test_train_split/training_set/')
trainSplitFiles = [coffDir + 'single_test_train_split/training_set/' + f for f in trainSplitFiles]

for f in trainSplitFiles:
    if 'embeddingKeys' not in f:
        with open(f, 'r') as inFH:
            for line in inFH:
                line = line.split(',')
                if line[0] != 'glycan':
                    line[0] = re.sub('"', '', line[0])
                    coffGlys.append(line[0])
            
coffGlys = list(set(coffGlys))


missingGlys = [g for g in coffGlys if g not in cfgGly2cfgID.keys()]
print(len(missingGlys))
missingGlys



# G-ol-Sp8
## Not in original RFU data, not sure why
cfgGly2cfgID['G-ol-Sp8'] = None 

# 'GlcNAcb1-4-MDPLys' --> 'GlcNAcb'
## Only GlcNAc monosaccharide without another defined linker
cfgGly2cfgID['GlcNAcb']
cfgGly2cfgID['GlcNAcb-Sp0']
getEmb(16)

cfgGly2cfgID['GlcNAcb1-4-MDPLys'] = cfgGly2cfgID['GlcNAcb']
getEmb(cfgGly2cfgID['GlcNAcb1-4-MDPLys'])


# Galb1-4GlcNAcb1-6(Galb1-4GlcNAcb1-2)Mana1-6(GlcNAcb1-4Galb1-4GlcNAcb1-4(Galb1-4GlcNAcb1-2)Mana1-3)Manb1-4GlcNAcb1-4(Fuca1-6)GlcNAc-Sp21
tmp = [k for k,v in cfgGly2cfgID.items() if 509 == v][0]
missingGlys[0]
tmp
## Original list of CFG glycans was missing a branch parenthesis, no embedding available
cfgGly2cfgID[missingGlys[0]] = None


missingGlys = [g for g in coffGlys if g not in cfgGly2cfgID.keys()]
print(len(missingGlys))
missingGlys



# Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-6(Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb-Sp21
tmp = [k for k,v in cfgGly2cfgID.items() if 57 == v][0]
tmp
missingGlys[0]
## Extra hyphen from iupac name original string, same glycan
cfgGly2cfgID[missingGlys[0]] = cfgGly2cfgID[tmp]
getEmb(57)


missingGlys = [g for g in coffGlys if g not in cfgGly2cfgID.keys()]
print(len(missingGlys))
missingGlys


# Gala1-3Fuca1-2Galb1-4GlcNAcb1-6(Gala1-3Fuca1-2Galb1-4GlcNAcb1-3)GalNAc-Sp14
## Matches "Gala1-3(Fuca1-2)Galb1-4GlcNAcb1-6(Gala1-3(Fuca1-2)Galb1-4GlcNAcb1-3)GalNAc-Sp14" from CFG glycans exactly besides the branching
tmp = [k for k,v in cfgGly2cfgID.items() if 451 == v][0]
missingGlys[0]
tmp
re.sub(r'[(|)]', '', missingGlys[0]) == re.sub(r'[(|)]', '', tmp)

tmp in coffGlys # Corresponding cfg glycan not otherwise in coff et al glycans
cfgGly2cfgID[missingGlys[0]] = cfgGly2cfgID[tmp]
getEmb(cfgGly2cfgID[missingGlys[0]])


missingGlys = [g for g in coffGlys if g not in cfgGly2cfgID.keys()]
print(len(missingGlys))
### All Coff et al glycans accounted for!

##############################

trainSplitFiles

trainSplitOut = [re.sub(r'DATA', 'embeddingKeys', f) for f in  trainSplitFiles]
trainSplitOut = [re.sub(r'.csv', r'.tsv', f) for f in  trainSplitOut]


for fIN, fOUT in zip(trainSplitFiles, trainSplitOut):
    with open(fIN, 'r') as fhIn:
        with open(fOUT, 'w') as fhOut:
            




































