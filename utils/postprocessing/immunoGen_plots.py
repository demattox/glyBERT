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
import networkx as nx

import cbk.sugarTrees.bin.sugarTrees as SugarTrees
import cbk.sugarTrees.bin.getStats as SugarStats

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)




homeDir = '/Users/dmattox/cbk/sugarTrees/'
# homeDir = '/Users/f002tsy/sugarTrees/'

# Load data
with open(homeDir + 'immunoFiles/allimmu_gen.pkl', 'rb') as inFH:
    allemb3 = pickle.load(inFH)
   
with open(homeDir + 'immunoFiles/allimmu_gen_label.pkl', 'rb') as inFH:
    allemblabel3 = pickle.load(inFH)
    
with open(homeDir + 'immunoFiles/npp.pkl', 'rb') as inFH: # Encoding of the generated glycans, overlap with allemblabel3
    npp = pickle.load(inFH)

with open(homeDir + 'immunoFiles/pls.pkl', 'rb') as inFH:
    pls = pickle.load(inFH)
pls = pls[::-1] # list of sizes of blocks of glycans generated along the same path, originally in reverse order starting at the backend of the dataframe

with open(homeDir + 'immunoFiles/path_encoding.pkl', 'rb') as inFH:
    encs = pickle.load(inFH)
encs = encs[::-1]
encs = [sublst[::-1] for sublst in encs] # lists of lists of the monosaccharides comprising the glycans at each step of each generative path

with open(homeDir + 'immunoFiles/id-path.pkl', 'rb') as inFH:
    pathIDs = pickle.load(inFH)

with open(homeDir + 'immunoFiles/probofemb.pkl', 'rb') as inFH:
    immProb = pickle.load(inFH) # Probabilities of being immunogenic/not immunogenic at each step of each generative path
# 1st col : probability of being not immunogenic
# 2nd col : probability of being immunogenic
immProb[:,0] + immProb[:,1]
# Same order as encs
immProb  = immProb[::-1,1] # Reversed to match embedding order, using only the immunogenic column

monoDict = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/commonMonoDict.pickle', 'rb'))
revMonoDict = {v:k for k,v in monoDict.items()}

with open(homeDir + 'pickles/matchedTrees_andLinks_andAno_immuno.pickle', 'rb') as inFH:
    matchedTrees =  pickle.load(inFH)


gDict = {g.sbID : g for g in SugarStats.commonGlys}


immDict = {0 : 'Non-immunogenic',
           1 : 'Immunogenic',
           2 : 'Generated'}
           # 3 : 'Non-imm. \n -- Generative Start'}


# generative paths
len(pls) # 20 paths in total

genPaths = []

for i, pathLen in enumerate(pls): # index paths from 1-20 and rep for path length to match glycan embeddings
    genPaths.extend([i+1] * pathLen)
    
print(len(genPaths) == npp.shape[0])

genPaths = [0] * sum(np.array(allemblabel3) != 2) + genPaths # Add 0 to all other existing points

print(len(genPaths) == len(allemblabel3))

# First glycan in each generative path is the starting, exisiting glycan, set label as such
prev = 1
for i in range(1,len(allemblabel3)):
    if allemblabel3[i-1] == 2 and genPaths[i] != prev:
        prev = genPaths[i]
        allemblabel3[i-1] = 0

allemblabel3[-1] = 0


# Full immProbs - get list of length matching full embedding array
allImmProbs = [0] * (allemb3.shape[0] - npp.shape[0]) # rep 0 for all exisiting glycan embeddings not used in paths
allImmProbs += list(immProb)



random.seed(27)

reducer = umap.UMAP().fit(allemb3) # fit UMAP to known & generated glycans
# reducer = umap.UMAP().fit(allemb3[np.array(allemblabel3) != 2,:]) # fit UMAP to known glycans
uu = reducer.embedding_
# uuGen = reducer.transform(allemb3[np.array(allemblabel3) == 2,:]) # Fit generated glycans to UMAP embedding built around known glycans

# uuGenp = pd.DataFrame(uuGen, columns = ['UMAP1', 'UMAP2'])
# uuGenp['Glycan Type'] = [immDict[lab] for lab in allemblabel3 if lab == 2]




uup = pd.DataFrame(uu, columns = ['UMAP1', 'UMAP2'])
# uup['Glycan Type'] = [immDict[lab] for lab in allemblabel3 if lab != 2]
uup['Glycan Type'] = [immDict[lab] for lab in allemblabel3]


# uup = pd.concat([uup, uuGenp])

p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Glycan Type') + p9.geom_point(alpha = 0.5, size = 3, color = '#ffffff') + \
    p9.scale_fill_manual(values = ['magenta','#a84c3e', '#3e6ca8'])




uup['Generative Path'] = [str(i) for i in genPaths]
# uup['Generative Path'] = genPaths

pathCols = [('#%06X' % random.randint(0, 0xFFFFFF)) for i in range(len(pls))]


immPathPlt = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Glycan Type', color = 'Generative Path', group = 'Generative Path') + \
    p9.geom_point(alpha = 0.5, size = 3, color = '#ffffff') + \
        p9.scale_fill_manual(values = ['magenta','#a84c3e', '#3e6ca8']) + \
            p9.geom_path() + \
                p9.scale_color_manual(values=['#FFFFFF00'] + pathCols, guide=False) + \
                    p9.theme_minimal()

immPathPlt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/immGenPaths.pdf', dpi = 600)

# add in probability of immunogenicity
uup['Immuno Probability'] = allImmProbs

# add in monosaccharide info
encsFlat = [subLst for lst in encs for subLst in lst]
len(encsFlat)

allEncs = [[] for i in range(allemb3.shape[0] - npp.shape[0])] # rep 0 for all exisiting glycan embeddings not used in paths
allEncs += encsFlat

uup['Monosaccs'] = allEncs


# Try and match paths to starting glycans
for sbGly in pathIDs.keys():
    if pathIDs[sbGly]:
        print(sbGly)
        for p in pathIDs[sbGly]:
            print('\t', p)
        print()

for i in range(1,21):
    print(uup[uup['Generative Path'] == str(i)])
    print('\n\n')

path2glyID = {'SBID6185': ['18','19', '20'],
 'SBID7112': ['16', '17'],
 'SBID10624': ['15'],
 'SBID12056': ['12', '13', '14'],
 'SBID13034': ['7', '8', '9', '10', '11'],
 'SBID17909': ['4', '5', '6'],
 'SBID4739': ['1', '2', '3']}
# Could've done this ^^ programatically but didn't realize it was ordered so nicely until I started doing it manually

len(pathIDs)

for p in pathIDs.keys():
    print(gDict[p].iupac)
    print(p, [d for d in matchedTrees if p in d['No']], 'paths:', path2glyID.get(p))

dLst = [d for p in pathIDs.keys() for d in matchedTrees if p in d['No']]
dLstNonGlys = set([g for d in dLst for g in d['No']])
dLstImmGlys = set([g for d in dLst for g in d['Yes']])

len(dLstNonGlys)
len(dLstImmGlys)



startingGlys = []
for p in genPaths:
    if p == 0:
        startingGlys.append('')
    else:
        for k,v in path2glyID.items():
            if str(p) in v:
                sGly = k
        startingGlys.append(sGly)

uup['Starting Glycan'] = startingGlys

# Update monosaccs with actual monosaccharide names
for i in range(len(allEncs)):
    if allEncs[i]:
        newEnc = [revMonoDict[m-2] for m in allEncs[i]]
        allEncs[i] = newEnc[::-1]

uup['Monosaccs'] = allEncs

sGly2iupac = {}
for k in path2glyID.keys():
    print(k)
    s = [sug for sug in SugarStats.commonGlys if sug.sbID == k][0]
    s.print()
    print('\n\n---------------\n')
    sGly2iupac[k] = s.iupac


for i in range(1,21):
    print(uup[uup['Generative Path'] == str(i)])
    print('\n_________________\n')

# Find IUPAC strings for each glycan
## Not built to handle branched glycans as all test cases here are linear glycans
iupacs = []
for monoLst, sGly in zip(allEncs, startingGlys):
    if sGly != '':
        sGly = sGly2iupac[sGly]
        decomposedGly = re.split(r'(\([a,b][0-9]-[0-9]\))', sGly)
        for mInd,gInd in zip(range(len(monoLst)), range(0, len(decomposedGly), 2)):
            decomposedGly[gInd] = monoLst[mInd]
        iupacs.append(''.join(decomposedGly))
    else:
        iupacs.append('')

uup['IUPAC'] = iupacs

# Check which generated glycans exist as labelled glycans
## Not built to handle branched glycans as all test cases here are linear glycans
knownGlycans = []
for gly,gType in zip(iupacs,allemblabel3):
    if gly != '':
        # print(gly, '\n\t', immDict[gType], '\n')
        s = [sug for sug in SugarStats.commonGlys if sug.iupac == gly]
        if s:
            # s[0].print()
            if len(s)> 1:
                print('Warning: Multiple matches, breaking loop...')
                break
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
            # print('NOVEL')
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
                
            # print('----------------------')
    else:
        knownGlycans.append('')

knownGlycans = [g if g != 'No' else 'Non-immunogenic' for g in knownGlycans]
knownGlycans = [g if g != 'Yes' else 'Immunogenic' for g in knownGlycans]


uup['Existing glycan'] = knownGlycans


uup[uup['Starting Glycan'] == 'SBID4739']
matches = [d for d in matchedTrees if 'SBID4739' in d['No']][0]
[gDict[g].print() for g in matches['Yes']]





with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/immGenUMAPall.pkl', 'wb') as outFH:
    pickle.dump(uup, outFH)

uup.to_csv(homeDir + 'immGenDF.tsv', sep="\t", index=False)

##############################################
# Draw network

for i in range(1,21):
    print(uup[uup['Generative Path'] == str(i)])
    print('\n_________________\n')

allGenGlys = list(set(uup[uup['IUPAC'] != '']['IUPAC']))
len(allGenGlys)

set(uup[uup['Generative Path'] == '19']['IUPAC'])

G = nx.DiGraph()
G.add_nodes_from(set(uup[uup['Generative Path'] == '19']['IUPAC']))
G.add_weighted_edges_from([('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'GalNAc(a2-3)Gal(b1-3)GalNAc', 1.0),
                           ('GalNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNAc(a2-3)Gal(b1-3)GalNAc', 1.0),
                           ('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNAc(a2-3)Gal(b1-3)GlcNAc', 1.0)])
nx.draw_networkx(G, node_color='r', edge_color='b',
                 pos = {'NeuNAc(a2-3)Gal(b1-3)GalNAc' : (0.852186,0),
                        'GalNAc(a2-3)Gal(b1-3)GalNAc' : (0.855664, 1),
                        'NeuNAc(a2-3)Gal(b1-3)GlcNAc' : (0.993761, -1)})

G = nx.DiGraph()
G.add_nodes_from([2, 3])








#####################

encs # monosacc code here -2 = monoDict key ; sublists reordered from imm to non-imm to match array rows (reverse from rows of array)
pls


i=19
print(encs[i-1])
print(uup[np.array(genPaths) == i])

p9.ggplot(uup[np.array(genPaths) == i]) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Glycan Type', color = 'Generative Path', group = 'Generative Path') + p9.geom_point(alpha = 0.5, size = 3, color = '#ffffff') + \
    p9.scale_fill_manual(values = ['magenta','#3e6ca8']) + \
        p9.geom_path()

# fig, ax = plt.subplots()
# ax.scatter(uup[np.array(genPaths) == i]['UMAP1'], uup[np.array(genPaths) == i]['UMAP2'])
# for j,pos in enumerate(zip(uup[np.array(genPaths) == i]['UMAP1'], uup[np.array(genPaths) == i]['UMAP2'])):
#     ax.annotate(str(j+1), xy = pos)
# plt.show()



changeDict = collections.defaultdict(list)
for p in range(1,21):
    glys = encs[p-1]
    print(p, '\n', glys)
    maxDelt = 0
    maxGlyDiff = ()
    for i in range(1,len(glys)):
        for m1,m2 in zip(glys[i],glys[i-1]):
            if m1 != m2:
                diff = (m1,m2)
        delt = np.sqrt(sum((uup[np.array(genPaths) == p].iloc[i-1, [0,1]] - uup[np.array(genPaths) == p].iloc[i, [0,1]])**2))
        distFromNonImmStart = [np.sqrt(sum((uup[np.array(genPaths) == p].iloc[-1, [0,1]] - uup[np.array(genPaths) == p].iloc[i, [0,1]])**2)),
                               np.sqrt(sum((uup[np.array(genPaths) == p].iloc[-1, [0,1]] - uup[np.array(genPaths) == p].iloc[i-1, [0,1]])**2))]
        direction = 1 if np.argmax(distFromNonImmStart) == 1 else -1
        print('\t\t',revMonoDict[diff[0]-2], '-->', revMonoDict[diff[1]-2], delt * direction)
        changeDict[diff].append(delt*direction)
        if delt > maxDelt:
            maxDelt = delt
            maxGlyDiff = diff
    print('\t',revMonoDict[maxGlyDiff[0]-2], '-->', revMonoDict[maxGlyDiff[1]-2], maxDelt * direction)
    print()

changeDict

with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/immunoFiles/glyChanges.csv', 'w') as outFH:
    outFH.write('old,new,distance\n')
    for monoPair, distLst in changeDict.items():
        for d in distLst:
            outFH.write(','.join([revMonoDict[monoPair[0] - 2], revMonoDict[monoPair[1] - 2], str(d)]) + '\n')









for k,v in changeDict.items():
    print(revMonoDict[k[0]-2], '-->', revMonoDict[k[1]-2])
    print(sorted(v))
    print()

meanChanges = [(k, np.mean(v)) for k,v in changeDict.items()]
meanChanges.sort(key = lambda x: x[1])
meanChanges

for k,score in meanChanges:
    print(revMonoDict[k[0]-2], '-->', revMonoDict[k[1]-2], score)

changeDictNew = collections.defaultdict(list)
for k,v in changeDict.items():
    changeDictNew[k[1]].extend(sorted(v))

meanChangesNew = [(k, np.mean(v)) for k,v in changeDictNew.items()]
meanChangesNew.sort(key = lambda x: x[1])
meanChangesNew


glyDiff = []
dists = []
for k,v in meanChanges:
    glyDiff.extend([revMonoDict[k[0]-2] + '>' + revMonoDict[k[1]-2]] * len(changeDict[k]))
    # glyDiff.extend([(revMonoDict[k[0]-2], revMonoDict[k[1]-2])] * len(changeDict[k]))
    dists.extend(sorted(changeDict[k]))

gdp = pd.DataFrame(list(zip(glyDiff, dists)), columns = ['Change', 'Distance'])
# gdp = pd.DataFrame(glyDiff, columns = ['Old', 'New'])
# gdp['Distance'] = dists

# gdpMean = 

p9.ggplot(gdp) + p9.aes(x = 'Change', y = 'Distance', fill = 'Distance') + \
    p9.geom_bar(stat = 'identity', position = 'dodge', color = 'black') + \
        p9.theme(axis_text_x = p9.element_text(angle = 45, vjust = 1, hjust=1, face = 'bold')) + \
            p9.scale_x_discrete(limits = [revMonoDict[g[0]-2]+'>'+revMonoDict[g[1]-2] for g,s in meanChanges]) + \
                p9.scale_fill_gradient(high = '#a84c3e', low = '#3e6ca8')












