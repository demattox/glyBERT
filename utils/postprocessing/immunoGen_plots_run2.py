#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 13:46:35 2021

@author: dmattox
"""

import pickle
import umap
import random
import collections
import re
import itertools

import plotnine as p9
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import networkx as nx

import cbk.sugarTrees.bin.sugarTrees as SugarTrees
import cbk.sugarTrees.bin.getStats as SugarStats

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)




homeDir = '/Users/dmattox/cbk/sugarTrees/'
# homeDir = '/Users/f002tsy/sugarTrees/'

# Load data

with open(homeDir + 'immunoFiles/run2/matchmono.pkl', 'rb') as inFH:
    matchmono = pickle.load(inFH)

with open(homeDir + 'immunoFiles/run2/notmatchmono.pkl', 'rb') as inFH:
    nomatchmono = pickle.load(inFH)

with open(homeDir + 'immunoFiles/id-path.pkl', 'rb') as inFH:
    r1_pathIDs = pickle.load(inFH)

with open(homeDir + 'pickles/matchedTrees_andLinks_andAno_immuno.pickle', 'rb') as inFH:
    matchedTrees =  pickle.load(inFH)


# with open(homeDir + 'immunoFiles/allimmu_gen.pkl', 'rb') as inFH:
#     allemb3 = pickle.load(inFH)
   
# with open(homeDir + 'immunoFiles/allimmu_gen_label.pkl', 'rb') as inFH:
#     allemblabel3 = pickle.load(inFH)
    
# with open(homeDir + 'immunoFiles/npp.pkl', 'rb') as inFH: # Encoding of the generated glycans, overlap with allemblabel3
#     npp = pickle.load(inFH)

# with open(homeDir + 'immunoFiles/pls.pkl', 'rb') as inFH:
#     pls = pickle.load(inFH)
# pls = pls[::-1] # list of sizes of blocks of glycans generated along the same path, originally in reverse order starting at the backend of the dataframe

# with open(homeDir + 'immunoFiles/path_encoding.pkl', 'rb') as inFH:
#     encs = pickle.load(inFH)
# encs = encs[::-1]
# encs = [sublst[::-1] for sublst in encs] # lists of lists of the monosaccharides comprising the glycans at each step of each generative path

# with open(homeDir + 'immunoFiles/id-path.pkl', 'rb') as inFH:
#     pathIDs = pickle.load(inFH)

# with open(homeDir + 'immunoFiles/probofemb.pkl', 'rb') as inFH:
#     immProb = pickle.load(inFH) # Probabilities of being immunogenic/not immunogenic at each step of each generative path
# # 1st col : probability of being not immunogenic
# # 2nd col : probability of being immunogenic
# immProb[:,0] + immProb[:,1]
# # Same order as encs
# immProb  = immProb[::-1,1] # Reversed to match embedding order, using only the immunogenic column

monoDict = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/commonMonoDict.pickle', 'rb'))
revMonoDict = {v:k for k,v in monoDict.items()}

with open(homeDir + 'pickles/matchedTrees_andLinks_andAno_immuno.pickle', 'rb') as inFH:
    matchedTrees =  pickle.load(inFH)


gDict = {g.sbID : g for g in SugarStats.commonGlys}


print(len(matchmono.keys()))

print(len(r1_pathIDs.keys()))


dropLst = []
for g in matchmono.keys():
    print(g)
    print('\t', gDict[g].iupac,'\n\t',len(matchmono[g]))
    if g in r1_pathIDs.keys():
        print('\tRun 1: ', len(r1_pathIDs[g]))
    structs = []
    structs = [s for s in matchedTrees if g in s['No']]
    if structs:
        print('Non-imm:', structs[0]['No'])
        print('Imm:', structs[0]['Yes'])
    else:
        print('Not in the matched trees,links, and anos dict\nRemove path from consideration')
        dropLst.append(g)
    print('\n\n')

for k in dropLst:
    del matchmono[k]

print(len(matchmono))



pIDs = [] # arbitrary numbering of generative paths spanning all starting points
startGlys = [] # Non-immunogenic starting structure for each path
immProb = [] # probability of immuonogenicity of each structure
mIDs = [] # lists of glycan IDs matching monoDict entries
isStart = [] # indicate the known non-immunogenic starting glycans
glyIUPACs = [] # Hold the iupac of the starting glycans

p=0 # initialize pID counter (first path gets counted at one)
for sGly in matchmono.keys():
    sGlyIUPAC = gDict[sGly].iupac
    for path in matchmono[sGly]:
        p+=1
        for i,g in enumerate(path): # For each structure along the generative path
            pIDs.append(p)
            startGlys.append(sGly)
            immProb.append(g[1][0][1])
            mIDs.append(g[0])
            if i == 0: # The first glycan in the path (list) is the starting point
                isStart.append(True)
                glyIUPACs.append(sGlyIUPAC)
            else:
                isStart.append(False)
                glyIUPACs.append('')
        


uup = pd.DataFrame({'Path' : pIDs,
                    'Start' : startGlys,
                    'Prob_imm' : immProb,
                    'isStart' : isStart,
                    'IUPAC' : glyIUPACs,
                    'MonoIDs' : mIDs})

# Update monosacc IDs with actual monosaccharide names
for i in range(len(mIDs)):
    if mIDs[i]:
        newEnc = [revMonoDict[m-2] for m in mIDs[i]]
        mIDs[i] = newEnc[::-1]

uup['Monosaccs'] = mIDs

# uup[['Path', 'IUPAC', 'Monosaccs']]


for i in range(1,10):
    print(uup[uup['Path'] == i])
    print('\n_________________\n')


################
# Double check for branched glycans, will need to be handled carefully when building IUPACS for generated glycans, or checking if they exist
uup[uup.apply(lambda row: '[' in row.IUPAC, axis = 1)]
uup[uup['Start'] == 'SBID10375']

# One starting glycan with a branch (SBID10375), with one generative path (63)
## Monosacc list order doesn't match IUPAC string order, swap 1st and 2nd elements in this list
swapLst = [[x.pop(1)]+x for x in uup[uup['Start'] == 'SBID10375']['Monosaccs'].to_list()]

for i,s in zip([239,240,241,242], swapLst):
    uup.at[i, 'Monosaccs'] = s


####
# Limit to glycans used in figure 4
####
print(uup.shape)
uup = uup.drop([239,240,241,242], axis=0)
print(uup.shape)

################



# Find IUPAC strings for each glycan
### Can only handle limited branching where the monosaccharide list is ordered to match the branch positions when split at each linkage

iupacs = []
for monoLst, sGly in zip(uup['Monosaccs'], uup['Start']):
    # print(sGly)
    sGly = gDict[sGly].iupac
    decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', sGly)
    for mInd,gInd in zip(range(len(monoLst)), range(0, len(decomposedGly), 2)):
        decomposedGly[gInd] = monoLst[mInd]
    iupacs.append(''.join(decomposedGly))
uup['IUPAC'] = iupacs


# Check which generated glycans exist as labelled glycans
## Not built to handle branched glycans as all test cases here are linear glycans
### Need to double check branched glycans in path 63 don't exist with the order of the branches flipped
knownGlycans = []
for gly,p in zip(uup['IUPAC'], uup['Path']):
    # print(gly, '\n\t', immDict[gType], '\n')
    if p == 63:
        print(gly)
    s = [sug for sug in SugarStats.commonGlys if sug.iupac == gly]
    if s:
        if len(s)> 1:
            print('Warning: Multiple matches, breaking loop...')
            break
        if p == 63:
            s[0].print()
        if s[0].immunogenic != 'Unknown':
            knownGlycans.append(s[0].immunogenic)
        else: # Check if unlabeled glycans are close to any known immunogenic glycans (same monos and links but different anomeric state)
            decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
            glyProds = []
            for m in decomposedGly:
                if re.match(r'\([a,b][0-9]-[0-9]\)', m):
                    out = [re.sub(r'\(b', '(a', m)]
                    out.append(re.sub(r'\(a', '(b', m))
                    glyProds.append(out)
                else:
                    glyProds.append([m])
            glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
            matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
            matchedLabs = []
            for s in matchedProds:
                if s:
                    matchedLabs.append(s[0].immunogenic)
            if 'Yes' in matchedLabs:
                knownGlycans.append('Unknown - Imm ano mismatch')
            else:
                knownGlycans.append('Unknown')
    else:
        if p == 63:
            print('NOVEL')
        # If novel, check if close to any known immunogenic glycans
        decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
        glyProds = []
        for m in decomposedGly:
            if re.match(r'\([a,b][0-9]-[0-9]\)', m):
                out = [re.sub(r'\(b', '(a', m)]
                out.append(re.sub(r'\(a', '(b', m))
                glyProds.append(out)
            else:
                glyProds.append([m])
        glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
        matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
        matchedLabs = []
        for s in matchedProds:
            if s:
                matchedLabs.append(s[0].immunogenic)
        if 'Yes' in matchedLabs:
            knownGlycans.append('Novel - Imm ano mismatch')
        else:
            knownGlycans.append('Novel')
    if p ==63:
        print('----------------------')


# ## Two of the branched glycans appear novel, check with flipped order in iupac
# # Fuc(a1-2)[Gal(a1-3)]Gal(b1-3)Fuc
# s = [sug for sug in SugarStats.commonGlys if sug.iupac == 'Gal(a1-3)[Fuc(a1-2)]Gal(b1-3)Fuc']
# print(s) # Doesn't appear in this order either, check anomeric isomers
# gly = 'Gal(a1-3)[Fuc(a1-2)]Gal(b1-3)Fuc'
# decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
# glyProds = []
# for m in decomposedGly:
#     if re.match(r'\([a,b][0-9]-[0-9]\)', m):
#         out = [re.sub(r'\(b', '(a', m)]
#         out.append(re.sub(r'\(a', '(b', m))
#         glyProds.append(out)
#     else:
#         glyProds.append([m])
# glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
# matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
# matchedLabs = []
# for s in matchedProds:
#     if s:
#         matchedLabs.append(s[0].immunogenic)
# print(matchedLabs)
# # Nope, definitely novel

# # Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-3)GlcNAc
# s = [sug for sug in SugarStats.commonGlys if sug.iupac == 'GalNAc(a1-3)[Fuc(a1-2)]Gal(b1-3)GlcNAc']
# print(s) # Not novel, known immunogenic!
# knownGlycans[-1] = s[0].immunogenic




knownGlycans = [g if g != 'No' else 'Non-immunogenic' for g in knownGlycans]
knownGlycans = [g if g != 'Yes' else 'Immunogenic' for g in knownGlycans]


uup['Existing glycan'] = knownGlycans


uup.to_csv(homeDir + 'immunoFiles/run2/matchmonoDF.tsv', sep="\t", index=False)






#################################################################
# Repeat with notmatchmono


len(nomatchmono.keys())

len(r1_pathIDs.keys())


dropLst = []
for g in nomatchmono.keys():
    print(g)
    print('\t', gDict[g].iupac,'\n\t',len(nomatchmono[g]))
    if g in r1_pathIDs.keys():
        print('\tRun 1: ', len(r1_pathIDs[g]))
    structs = []
    structs = [s for s in matchedTrees if g in s['No']]
    if structs:
        print('Non-imm:', structs[0]['No'])
        print('Imm:', structs[0]['Yes'])
    else:
        print('Not in the matched trees,links, and anos dict\nRemove path from consideration')
        dropLst.append(g)
    print('\n\n')

for k in dropLst:
    del nomatchmono[k]

print(len(nomatchmono))



pIDs = [] # arbitrary numbering of generative paths spanning all starting points
startGlys = [] # Non-immunogenic starting structure for each path
immProb = [] # probability of immuonogenicity of each structure
mIDs = [] # lists of glycan IDs matching monoDict entries
isStart = [] # indicate the known non-immunogenic starting glycans
glyIUPACs = [] # Hold the iupac of the starting glycans

p=0 # initialize pID counter (first path gets counted at one)
for sGly in nomatchmono.keys():
    sGlyIUPAC = gDict[sGly].iupac
    for path in nomatchmono[sGly]:
        p+=1
        for i,g in enumerate(path): # For each structure along the generative path
            pIDs.append(p)
            startGlys.append(sGly)
            immProb.append(g[1][0][1])
            mIDs.append(g[0])
            if i == 0: # The first glycan in the path (list) is the starting point
                isStart.append(True)
                glyIUPACs.append(sGlyIUPAC)
            else:
                isStart.append(False)
                glyIUPACs.append('')
        


uup = pd.DataFrame({'Path' : pIDs,
                    'Start' : startGlys,
                    'Prob_imm' : immProb,
                    'isStart' : isStart,
                    'IUPAC' : glyIUPACs,
                    'MonoIDs' : mIDs})

# Update monosacc IDs with actual monosaccharide names
for i in range(len(mIDs)):
    if mIDs[i]:
        newEnc = [revMonoDict[m-2] for m in mIDs[i]]
        mIDs[i] = newEnc[::-1]

uup['Monosaccs'] = mIDs

# uup[['Path', 'IUPAC', 'Monosaccs']]


for i in range(1,10):
    print(uup[uup['Path'] == i])
    print('\n_________________\n')


################
# Double check for branched glycans, will need to be handled carefully when building IUPACS for generated glycans, or checking if they exist
uup.index[uup.apply(lambda row: '[' in row.IUPAC, axis = 1)]
branchedGlys = uup.index[uup['Start'] == 'SBID10375']

# # One starting glycan with a branch (SBID10375), with one generative path (63)
# ## Monosacc list order doesn't match IUPAC string order, swap 1st and 2nd elements in this list
# swapLst = [[x.pop(1)]+x for x in uup[uup['Start'] == 'SBID10375']['Monosaccs'].to_list()]

# for i,s in zip([239,240,241,242], swapLst):
#     uup.at[i, 'Monosaccs'] = s

# So many paths and glys, Fuc branching isn't one of the original 8, leave out here
print(uup.shape)
uup = uup.drop(branchedGlys, axis=0)
print(uup.shape)
################



# Find IUPAC strings for each glycan
### Can only handle limited branching where the monosaccharide list is ordered to match the branch positions when split at each linkage

iupacs = []
for monoLst, sGly in zip(uup['Monosaccs'], uup['Start']):
    # print(sGly)
    sGly = gDict[sGly].iupac
    decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', sGly)
    for mInd,gInd in zip(range(len(monoLst)), range(0, len(decomposedGly), 2)):
        decomposedGly[gInd] = monoLst[mInd]
    iupacs.append(''.join(decomposedGly))
uup['IUPAC'] = iupacs


# Check which generated glycans exist as labelled glycans
## Not built to handle branched glycans as all test cases here are linear glycans
### Need to double check branched glycans in path 63 don't exist with the order of the branches flipped
knownGlycans = []
for gly,p in zip(uup['IUPAC'], uup['Path']):
    # print(gly, '\n\t', immDict[gType], '\n')
    if p == 63:
        print(gly)
    s = [sug for sug in SugarStats.commonGlys if sug.iupac == gly]
    if s:
        if len(s)> 1:
            print('Warning: Multiple matches, breaking loop...')
            break
        if p == 63:
            s[0].print()
        if s[0].immunogenic != 'Unknown':
            knownGlycans.append(s[0].immunogenic)
        else: # Check if unlabeled glycans are close to any known immunogenic glycans (same monos and links but different anomeric state)
            decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
            glyProds = []
            for m in decomposedGly:
                if re.match(r'\([a,b][0-9]-[0-9]\)', m):
                    out = [re.sub(r'\(b', '(a', m)]
                    out.append(re.sub(r'\(a', '(b', m))
                    glyProds.append(out)
                else:
                    glyProds.append([m])
            glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
            matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
            matchedLabs = []
            for s in matchedProds:
                if s:
                    matchedLabs.append(s[0].immunogenic)
            if 'Yes' in matchedLabs:
                knownGlycans.append('Unknown - Imm ano mismatch')
            else:
                knownGlycans.append('Unknown')
    else:
        if p == 63:
            print('NOVEL')
        # If novel, check if close to any known immunogenic glycans
        decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
        glyProds = []
        for m in decomposedGly:
            if re.match(r'\([a,b][0-9]-[0-9]\)', m):
                out = [re.sub(r'\(b', '(a', m)]
                out.append(re.sub(r'\(a', '(b', m))
                glyProds.append(out)
            else:
                glyProds.append([m])
        glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
        matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
        matchedLabs = []
        for s in matchedProds:
            if s:
                matchedLabs.append(s[0].immunogenic)
        if 'Yes' in matchedLabs:
            knownGlycans.append('Novel - Imm ano mismatch')
        else:
            knownGlycans.append('Novel')
    if p ==63:
        print('----------------------')


# ## Two of the branched glycans appear novel, check with flipped order in iupac
# # Fuc(a1-2)[Gal(a1-3)]Gal(b1-3)Fuc
# s = [sug for sug in SugarStats.commonGlys if sug.iupac == 'Gal(a1-3)[Fuc(a1-2)]Gal(b1-3)Fuc']
# print(s) # Doesn't appear in this order either, check anomeric isomers
# gly = 'Gal(a1-3)[Fuc(a1-2)]Gal(b1-3)Fuc'
# decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\)[\[,\]]?)', gly)
# glyProds = []
# for m in decomposedGly:
#     if re.match(r'\([a,b][0-9]-[0-9]\)', m):
#         out = [re.sub(r'\(b', '(a', m)]
#         out.append(re.sub(r'\(a', '(b', m))
#         glyProds.append(out)
#     else:
#         glyProds.append([m])
# glyProds = [''.join(g) for g in list(itertools.product(*glyProds))]
# matchedProds = [[sug for sug in SugarStats.commonGlys if sug.iupac == g] for g in glyProds]
# matchedLabs = []
# for s in matchedProds:
#     if s:
#         matchedLabs.append(s[0].immunogenic)
# print(matchedLabs)
# # Nope, definitely novel

# # Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-3)GlcNAc
# s = [sug for sug in SugarStats.commonGlys if sug.iupac == 'GalNAc(a1-3)[Fuc(a1-2)]Gal(b1-3)GlcNAc']
# print(s) # Not novel, known immunogenic!
# knownGlycans[-1] = s[0].immunogenic




knownGlycans = [g if g != 'No' else 'Non-immunogenic' for g in knownGlycans]
knownGlycans = [g if g != 'Yes' else 'Immunogenic' for g in knownGlycans]


uup['Existing glycan'] = knownGlycans


uup.to_csv(homeDir + 'immunoFiles/run2/nomatchmonoDF.tsv', sep="\t", index=False)
