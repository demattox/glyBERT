#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:37:08 2021

@author: dmattox
"""


import pickle
import collections
import Levenshtein

import numpy as np
import matplotlib.pyplot as plt

import cbk.sugarTrees.bin.sugarTrees as SugarTrees

homeDir = '/Users/dmattox/cbk/sugarTrees/'

# Load encodings containing common monosaccharides
with open(homeDir + 'pickles/commonGlys_encodings.pickle', 'rb') as inFH:
    commonGlys = pickle.load(inFH)



# Limit to glycans with immunogenicity labels
immGlys = [gly for gly in commonGlys if gly.immunogenic != 'Unknown']
len(immGlys)

nodeCnts = [len(gly.tree) for gly in immGlys]
print(collections.Counter(nodeCnts).most_common())

nodeCnts = list(set(nodeCnts))
# Look within sets of glycans with the same nuber of nodes for other glycans with identical tree structures


#######
# Get matched tree structures within immunogenicity-labelled glycans
#######

matchedTrees = collections.defaultdict(lambda: collections.defaultdict(list))
## Keys: list of dict keys for tree (unique for each tee structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in immGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = tuple(g.tree.keys())
        matchedTrees[treeKey][g.immunogenic].append(g.sbID)

print(len(immGlys))
print(len(matchedTrees))
 # In the 1127 glycans with immunogenicity labels, there are 364 unique tree structures
 
print(len([v for v in matchedTrees.values() if len(v) == 2]))
# 36 of those 364 tree structures have at least one glycan that is imm and at least one that is not

#########################
# Looking for glycans with similar structures and/or saccharide composition but different immunogenic labels
# nodeCnts = [len(gly.tree) for gly in immGlys]
# print(collections.Counter(nodeCnts).most_common())

# nonImmGlys = [g for g in immGlys if g.immunogenic == 'No']
# allImmGlys = [g for g in immGlys if g.immunogenic == 'Yes']


# [g for g in nonImmGlys if 'LDManHep' in g.iupac]
# [g for g in nonImmGlys if 'Kdo' in g.iupac]
# [g for g in nonImmGlys if 'Ara' in g.iupac]


# for nCnt in set(nodeCnts):
#     gList = []
#     for gly in immGlys:
#         if len(gly.tree) == nCnt:
#             gList.append(gly)
#     if len(set([g.immunogenic for g in gList])) > 1:
#         print(nCnt,collections.Counter(nodeCnts)[nCnt])

# gList = []
# for gly in immGlys:
#     if len(gly.tree) == 12:
#         gList.append(gly)
# [g.iupac for g in gList if g.immunogenic == 'Yes']

# [g.iupac for g in gList if g.immunogenic == 'No']


# [g for g in gList if g.immunogenic == 'Yes'][0].print()


# gList[-2].treePrint()
# gList[-2].iupac

# maxPair = (0,'')
# for g in gList:
#     if g.immunogenic == 'No':
#         print(g.sbID, '\n', g.iupac)
#         d = Levenshtein.ratio(g.iupac, gList[-2].iupac)
#         print('\t', d)
#         if d > maxPair[0]:
#             maxPair = (d, g)

# maxPair[1].treePrint()
# maxPair[1].iupac

# ###

# g1 = [g for g in immGlys if g.sbID == 'SBID15739'][0]
# g1.print()

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')
# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]

# ###

# g1 =[g for g in immGlys if g.sbID == 'SBID6861'][0]
# g1.print()

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')
# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]
# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# [g for g in nonImmGlys if g.iupac[-6:] == 'NeuNAc']
# [g for g in nonImmGlys if g.iupac == 'NeuNAc(a2-6)Glc(a1-6)[Glc(a1-4)Gal(b1-6)Glc(b1-3)]Gal(b1-3)[Glc(a1-6)]GalNAc']

# ###

# g1 =[g for g in immGlys if g.sbID == 'SBID409'][0]
# g1.print()

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-30,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')
# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]
# for gID in matchedTrees[k]['No']:
#     g = [g for g in nonImmGlys if g.sbID == gID][0]
#     g.print()
#     print('______')
# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# [g for g in nonImmGlys if g.iupac == 'GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[GalNAc(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc']




# for nG in nonImmGlys:
#     if '[Fuc(a1-3)]' in nG.iupac:
#         print(nG.iupac)
#         sub6 = [g for g in allImmGlys if g.iupac == nG.iupac.replace('[Fuc(a1-6)]', '')]
#         if sub6:
#             sub6[0].print()
#             print('\n')


# # Adding a1-3 fuc removes immunogenicity
# g1 = [g for g in immGlys if 'Fuc(a1-2)Gal(b1-4)GlcNAc' == g.iupac][0]
# g1.print()
# [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc'][0].print()


# for nG in nonImmGlys:
#     if '[Fuc(a1-3)]' in nG.iupac:
#         print(nG.iupac)
#         simLst = [(g, Levenshtein.ratio(g.iupac, nG.iupac)) for g in allImmGlys]
#         simLst.sort(key = lambda x: x[1])
#         ind = -1
#         while simLst[ind][1] > 0.9:
#             simLst[ind][0].print()
#             ind -= 1
#         print('---------------')

# g1 = [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc'][0]
# g1.print()
# [g for g in allImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc'][0].print()

# [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Glc'][0].print()
# [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)Glc'][0].print()
# [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc'][0].print()
# [g for g in nonImmGlys if g.iupac == 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)Glc'][0].print()


# for nG in nonImmGlys:
#         simLst = [(g, Levenshtein.ratio(g.iupac, nG.iupac)) for g in allImmGlys]
#         simLst.sort(key = lambda x: x[1])
#         ind = -1
#         if simLst[ind][1] > 0.9:
#             nG.print()
#             print('\n')
#         while simLst[ind][1] > 0.9:
#             simLst[ind][0].print()
#             ind -= 1
#         if simLst[-1][1] > 0.9:
#             print('------------------------------------------------------------\n------------------------------------------------------------')
#             null = input('press enter to continue')

# # Check fucose containing immunogenic glycans

# immFuc = []
# for g in immGlys:
#     if g.immunogenic == 'Yes' and 'Fuc' in g.iupac:
#         immFuc.append(g)

# len(immFuc)

# iGLst = []
# for iG in immFuc:
#     if '[Fuc(a1-3)]' in iG.iupac:
#         print(iG.iupac)
#         sub6 = [g for g in nonImmGlys if g.iupac == iG.iupac.replace('[Fuc(a1-3)]', '')]
#         if sub6:
#             sub6[0].print()
#             print('\n')
#             iGLst.append(iG)

# iGLst[0].print()
# [g for g in nonImmGlys if g.iupac == iGLst[0].iupac.replace('[Fuc(a1-3)]', '')][0].print()

# baseCnts = collections.Counter([n.base for n in iGLst[0].tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]
# compMatch[0].print()


# iGLst[1].print()
# [g for g in nonImmGlys if g.iupac == iGLst[1].iupac.replace('[Fuc(a1-3)]', '')][0].print()
# [g for g in nonImmGlys if g.iupac == iGLst[1].iupac.replace('[Fuc(a1-3)]', '[Fuc(a1-4)]')][0].print()


# baseCnts = collections.Counter([n.base for n in iGLst[1].tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]


# iGLst[2].print()
# [g for g in nonImmGlys if g.iupac == iGLst[2].iupac.replace('[Fuc(a1-3)]', '')][0].print()


# baseCnts = collections.Counter([n.base for n in iGLst[2].tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# ###
# g1 = [g for g in immFuc if 'Man(a1-6)[Xyl(b1-2)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-3)]GlcNAc' in g.iupac][0]
# g1.print()
# g1.species

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')

# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# for gID in matchedTrees[k]['No']:
#     g = [g for g in nonImmGlys if g.sbID == gID][0]
#     g.print()
#     print('______')

# [g for g in nonImmGlys if g.iupac == g1.iupac.replace('[Xyl(b1-2)]', '')][0].print()

# ###

# g1 = [g for g in immFuc if 'Gal(b1-4)[Fuc(a1-3)][Glc(a1-6)]GlcNAc(b1-3)Gal' in g.iupac][0]
# g1.print()
# g1.treePrint()

# g2 = [g for g in immFuc if 'Gal(b1-4)[Fuc(a1-3)][Gal(a1-6)]GlcNAc(b1-3)Gal' in g.iupac][0]

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')

# partials = [p for p in simLst if 'Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal' in p[0].iupac]
# for p in partials:
#     print(p)
#     p[0].print()
#     print('________________')

# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]
# compMatch[0].print()

# baseCnts = collections.Counter([n.base for n in g2.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]

# ###

# g1 = [g for g in immFuc if 'NeuNAc(a2-3)[Fuc(a1-2)Gal(b1-3)GalNAc(b1-4)]Gal(b1-4)Glc' in g.iupac][0]
# g1.print()

# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]
# for gID in matchedTrees[k]['Yes']:
#     g = [g for g in immGlys if g.sbID == gID][0]
#     g.print()
#     print('______')

# g2 = [g for g in immGlys if 'NeuNAc(a2-3)[NeuNAc(a2-3)Gal(b1-3)GalNAc(b1-4)]Gal(b1-4)Glc' in g.iupac][0]

# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')

# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]
# compMatch[0].print()

# baseCnts = collections.Counter([n.base for n in g2.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]

# ###

# g1 = [g for g in immFuc if 'NeuNAc(a2-6)Gal(b1-4)GlcNAc(b1-3)Fuc' in g.iupac][0]
# g1.print()

# k = [k for k in matchedTrees.keys() if g1.sbID in matchedTrees[k]['Yes']][0]
# matchedTrees[k]
# for gID in matchedTrees[k]['No']:
#     g = [g for g in immGlys if g.sbID == gID][0]
#     g.print()
#     print('______')
    
# simLst = [(g, Levenshtein.ratio(g.iupac, g1.iupac)) for g in nonImmGlys]
# simLst.sort(key = lambda x: x[1])
# for i in range(-10,-1,1):
#     print(simLst[i])
#     simLst[i][0].print()
#     print('________________')

# baseCnts = collections.Counter([n.base for n in g1.tree.values() if n.base not in ['START', 'END']])
# compMatch = [g for g in nonImmGlys if baseCnts == collections.Counter([n.base for n in g.tree.values() if n.base not in ['START', 'END']])]
# compMatch[0].print()

# [g for g in nonImmGlys if g.iupac == 'NeuNAc(a2-6)Gal(b1-4)GlcNAc'][0].print()

# [g for g in immGlys if g.iupac == 'NeuNAc(a2-3)Gal(b1-4)GlcNAc(b1-3)Fuc'][0].print()

# [g for g in immGlys if g.iupac == 'NeuNAc(a2-3)Gal(b1-4)GlcNAc'][0].print()


# g1 = [g for g in commonGlys if g.iupac == 'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc'][0]
# g1.print()
# g1.treePrint()
# for k,n in g1.tree.items():
#     print(k,n, n.subway)



# g1 = [g for g in commonGlys if g.iupac == 'NeuNAc(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc'][0]
# g1.print()
# g1.treePrint()
# for k,n in g1.tree.items():
#     print(k,n, n.subway, n.parentLink)


# g2 = [g for g in commonGlys if g.iupac == 'NeuNAc(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[NeuNAc(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc'][0]
# g2.print()
# g2.treePrint()
# for k,n in g2.tree.items():
#     print(k,n, n.subway, n.parentLink)


# ###

# nodeCnts = list(set(nodeCnts))
#########################

matchedTrees = [v for v in matchedTrees.values() if len(v) == 2]
matchedTrees

[g for g in immGlys if g.sbID == 'SBID409'][0].treePrint()
[g for g in immGlys if g.sbID == 'SBID409'][0].iupac
# [g for g in immGlys if g.sbID == 'SBID8965'][0].treePrint()
# [g for g in immGlys if g.sbID == 'SBID8965'][0].iupac

[g for g in immGlys if g.sbID == 'SBID4237'][0].treePrint()
[g for g in immGlys if g.sbID == 'SBID4237'][0].iupac


with open(homeDir + 'pickles/matchedTrees_immuno.pickle', 'wb') as outFH:
    pickle.dump(matchedTrees, outFH, protocol=pickle.HIGHEST_PROTOCOL)
    
    


#######
# Get matched tree structures (including C links) within immunogenicity-labelled glycans
#######

matchedLinkTrees = collections.defaultdict(lambda: collections.defaultdict(list))
## Keys: list of dict keys for tree (unique for each tree structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in immGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = list(g.tree.keys())
        treeKey.extend([tuple(v.Clinks.keys()) for v in g.tree.values()]) # Add tuples containing linkage info to the tree structure key
        matchedLinkTrees[tuple(treeKey)][g.immunogenic].append(g.sbID)
        
        if g.sbID == 'SBID7236':
            print(g.sbID, tuple(treeKey))
        if g.sbID == 'SBID1174':
            print(g.sbID, tuple(treeKey))


print(len(immGlys))
print(len(matchedLinkTrees))
 # In the 1127 glycans with immunogenicity labels, there are 801 unique tree structures
 
print(len([v for v in matchedLinkTrees.values() if len(v) == 2]))
# 31 of those 806 tree structures have at least one glycan that is imm and at least one that is not


matchedLinkTrees = [v for v in matchedLinkTrees.values() if len(v) == 2]

with open(homeDir + 'pickles/matchedTrees_andLinks_immuno.pickle', 'wb') as outFH:
    pickle.dump(matchedLinkTrees, outFH, protocol=pickle.HIGHEST_PROTOCOL)


nEx = [g for g in immGlys if g.sbID == 'SBID2031'][0]
pEx = [g for g in immGlys if g.sbID == 'SBID16265'][0]
nEx.immunogenic
pEx.immunogenic

nEx.treePrint()
[v.Clinks for v in nEx.tree.values()]
nEx.iupac

pEx.treePrint()
[v.Clinks for v in pEx.tree.values()]
pEx.iupac



#######
# Get matched tree structures (including anomeric state) within immunogenicity-labelled glycans
#######

matchedLinkAnoTrees = collections.defaultdict(lambda: collections.defaultdict(list))
## Keys: list of dict keys for tree (unique for each tree structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in immGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = list(g.tree.keys())
        treeKey.extend([tuple(v.Clinks.keys()) for v in g.tree.values()]) # Add tuples containing linkage info to the tree structure key
        treeKey.append(''.join([node.anomeric for node in g.tree.values()])) # Add string for unique anomeric conformations to the tree structure key
        matchedLinkAnoTrees[tuple(treeKey)][g.immunogenic].append(g.sbID)


print(len(immGlys))
print(len(matchedLinkAnoTrees))
 # In the 1127 glycans with immunogenicity labels, there are 898 unique tree structures
 
print(len([v for v in matchedLinkAnoTrees.values() if len(v) == 2]))
# 20 of those 898 tree structures have at least one glycan that is imm and at least one that is not


[v for v in matchedLinkAnoTrees.values() if len(v) == 2]


nEx = [g for g in immGlys if g.sbID == 'SBID6680'][0]
pEx = [g for g in immGlys if g.sbID == 'SBID8153'][0]
nEx.immunogenic
pEx.immunogenic

nEx.treePrint()
[v.Clinks for v in nEx.tree.values()]
nEx.iupac

pEx.treePrint()
[v.Clinks for v in pEx.tree.values()]
pEx.iupac


matchedLinkAnoTrees = [v for v in matchedLinkAnoTrees.values() if len(v) == 2]


with open(homeDir + 'pickles/matchedTrees_andLinks_andAno_immuno.pickle', 'wb') as outFH:
    pickle.dump(matchedLinkAnoTrees, outFH, protocol=pickle.HIGHEST_PROTOCOL)











#############
# Get all matched tree structures (ignore C links and anomeric conf.) from the entire set of 16,048 glycans
#############

nodeCnts = [len(gly.tree) for gly in commonGlys]
collections.Counter(nodeCnts).most_common()

nodeCnts = list(set(nodeCnts))

matchedTrees = collections.defaultdict(list)
## Keys: list of dict keys for tree (unique for each tee structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in commonGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = tuple(g.tree.keys())
        matchedTrees[treeKey].append(g.sbID)

print(len(commonGlys))
print(len(matchedTrees))
print(sum([1 for v in matchedTrees.values() if len(v) >= 2]))
 # In the 16048 glycans, there are 2112 unique tree structures, 970 of which have more than one glycan
 
treeMembers = [len(v) for v in matchedTrees.values()]
treeSizes = []
for treeStruct, idLst in matchedTrees.items():
    glyEx = [g for g in commonGlys if g.sbID == idLst[0]][0]
    treeSizes.append(sum([1 for node in glyEx.tree.values() if node not in ['START', 'END']]))


plt.scatter(treeSizes, np.log10(treeMembers), alpha = 0.5)
plt.xlabel('# of non-start/stop nodes')
plt.ylabel('log10(# of glycans with same structure)')


with open(homeDir + 'pickles/matchedTrees_all.pickle', 'wb') as outFH:
    pickle.dump(matchedTrees, outFH, protocol=pickle.HIGHEST_PROTOCOL)



#############
# Get all matching tree structures and glycosidic linkages (ignore ano) from all glycans
#############

nodeCnts = [len(gly.tree) for gly in commonGlys]
collections.Counter(nodeCnts).most_common()

nodeCnts = list(set(nodeCnts))

matchedLinkTrees_all = collections.defaultdict(list)
## Keys: list of dict keys for tree (unique for each tee structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in commonGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = list(g.tree.keys())
        treeKey.extend([tuple(v.Clinks.keys()) for v in g.tree.values()]) # Add tuples containing linkage info to the tree structure key
        matchedLinkTrees_all[tuple(treeKey)].append(g.sbID)
        
        if g.sbID == 'SBID7236':
            print(g.sbID, tuple(treeKey))
        if g.sbID == 'SBID1174':
            print(g.sbID, tuple(treeKey))

print(len(commonGlys))
print(len(matchedLinkTrees_all))
print(sum([1 for v in matchedLinkTrees_all.values() if len(v) >= 2]))
 # In the 16048 glycans, there are 6349 unique tree structures, 2107 of which have more than one glycan
 
treeMembers = [len(v) for v in matchedLinkTrees_all.values()]
treeSizes = []
for treeStruct, idLst in matchedLinkTrees_all.items():
    glyEx = [g for g in commonGlys if g.sbID == idLst[0]][0]
    treeSizes.append(sum([1 for node in glyEx.tree.values() if node not in ['START', 'END']]))

plt.scatter(treeSizes, np.log10(treeMembers), alpha = 0.5)
plt.xlabel('# of non-start/stop nodes')
plt.ylabel('log10(# of glycans with same structure)')


with open(homeDir + 'pickles/matchedTrees_andLinks_all.pickle', 'wb') as outFH:
    pickle.dump(matchedLinkTrees_all, outFH, protocol=pickle.HIGHEST_PROTOCOL)




with open(homeDir + 'pickles/matchedTrees_andLinks_all.pickle', 'rb') as inFH:
    tmp = pickle.load(inFH)
len(tmp.keys())





#############
# Get all matched tree structures, C links and anomeric conf. from the entire set of 16,048 glycans
#############

nodeCnts = [len(gly.tree) for gly in commonGlys]
collections.Counter(nodeCnts).most_common()

nodeCnts = list(set(nodeCnts))

matchedLinkAnoTrees_all = collections.defaultdict(list)
## Keys: list of dict keys for tree (unique for each tee structure)
## values: dict with two entries (imm and not imm), each entry is a dict([]) of SBIDs

for l in nodeCnts:
    glys = [gly for gly in commonGlys if len(gly.tree) == l]
    for g in glys:
        treeKey = list(g.tree.keys())
        treeKey.extend([tuple(v.Clinks.keys()) for v in g.tree.values()]) # Add tuples containing linkage info to the tree structure key
        treeKey.append(''.join([node.anomeric for node in g.tree.values()])) # Add string for unique anomeric conformations to the tree structure key
        matchedLinkAnoTrees_all[tuple(treeKey)].append(g.sbID)

print(len(commonGlys))
print(len(matchedLinkAnoTrees_all))
print(sum([1 for v in matchedLinkAnoTrees_all.values() if len(v) >= 2]))
 # In the 16048 glycans, there are 2112 unique tree structures, 970 of which have more than one glycan
 
treeMembers = [len(v) for v in matchedLinkAnoTrees_all.values()]
treeSizes = []
for treeStruct, idLst in matchedLinkAnoTrees_all.items():
    glyEx = [g for g in commonGlys if g.sbID == idLst[0]][0]
    treeSizes.append(sum([1 for node in glyEx.tree.values() if node not in ['START', 'END']]))


plt.scatter(treeSizes, np.log10(treeMembers), alpha = 0.5)
plt.xlabel('# of non-start/stop nodes')
plt.ylabel('log10(# of glycans with same structure)')


with open(homeDir + 'pickles/matchedTrees_andLinks_andAno_all.pickle', 'wb') as outFH:
    pickle.dump(matchedLinkAnoTrees_all, outFH, protocol=pickle.HIGHEST_PROTOCOL)
















