#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 13:11:46 2020

@author: dmattox
"""

import pickle

import cbk.sugarTrees.bin.sugarTrees as SugarTrees
import cbk.sugarTrees.bin.getStats as SugarStats

homeDir = '/Users/dmattox/cbk/sugarTrees/'

# sbFile = homeDir+'data/SugarBase_v2_2020_11_19.csv'

# SugarStats.commonGlys[0].print()
# SugarStats.commonGlys[0].treePrint()
# SugarStats.commonGlys[0].taxonomy
# SugarStats.commonGlys[0].encoding

SugarStats.maxB
SugarStats.maxC

set(SugarStats.ano)


 
# Check immunogenicity info
imm = set()
for sugar in SugarStats.commonGlys:
    imm.update([sugar.immunogenic.strip()])
print(imm)

immDict = {'Unknown': -1, 'No': 0, 'Yes': 1}

# check links
lis = set()
for sugar in SugarStats.commonGlys:
    lis.update([sugar.link.strip()])
print(lis)

linkDict = {'None': -1, 'Free': 0, 'N': 1, 'O': 2}



# colNames = ['sugar'] + ['anomeric'] + ['C_links'] + ['B_lines'] + ['parLink', 'sDepth', 'sInd']

  
#####
## Updated encoding (7 columns)
# Get dict to pickle with glycans with only common monosaccharides
glycan_encodings = {}

topColVals = [0]*7

for sugar in SugarStats.commonGlys:
    
    sugar.buildEncoding(SugarStats.maxB, SugarStats.maxC, SugarStats.monoDict)
    
    for i in range(len(topColVals)):
        if max(sugar.encoding[:,i]) > topColVals[i]:
            topColVals[i] = max(sugar.encoding[:,i])
            
        
    li = linkDict[sugar.link] # Hold link code
    if sugar.taxonomy:
        sp = []
        for i,r in enumerate(sugar.taxonomy):
            if r:
                sp.append([SugarStats.taxDicts[i][t] for t in r])
            else:
                sp.append([-1])
    else:
        sp = [[-1]] * 6
    # [skdict[s] for s in sugar.species] if sugar.taxonomy else [-1]*6 # Hold list of species codes
    im = immDict[sugar.immunogenic]
    
    glycan_encodings[sugar.sbID] = [sugar.encoding, li, sp, im]


print(len(glycan_encodings))

print(topColVals)

sugar.print()
sugar.encoding

glycan_encodings['SBID124']



# np.array_equal(sugar.encoding, glycan_encodings[sugar.sbID][0])





# encodeFile = homeDir+'pickles/glycan_encodings_6taxa5kings.pickle'

encodeFile = homeDir+'pickles/glycan_common_16048_may22.pickle'

with open(encodeFile, 'wb') as outFH:
    pickle.dump(glycan_encodings, outFH, protocol=pickle.HIGHEST_PROTOCOL)


# Pickled file read in with code below
# with open(encodeFile, 'rb') as inFH:
#     glycan_encodings = pickle.load(inFH)

# Pickle monosaccharide dictionary

mDictFile = homeDir + 'pickles/commonMonoDict.pickle'

with open(mDictFile, 'wb') as outFH:
    pickle.dump(SugarStats.monoDict, outFH, protocol=pickle.HIGHEST_PROTOCOL)


# Pickle reversed taxonomic id dictionaries
rev_taxDicts = []
for d in SugarStats.taxDicts:
    rev_taxDicts.append({v:k for k,v in d.items()})

taxid_num_dictFile = homeDir + 'pickles/taxid_num_dicts.pickle'
with open(taxid_num_dictFile, 'wb') as outFH:
    pickle.dump(rev_taxDicts, outFH, protocol=pickle.HIGHEST_PROTOCOL)
    
with open(homeDir + 'pickles/commonGlys_encodings.pickle', 'wb') as outFH:
    pickle.dump(SugarStats.commonGlys, outFH, protocol=pickle.HIGHEST_PROTOCOL)
    

