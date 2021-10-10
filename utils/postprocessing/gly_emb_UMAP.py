#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 13:32:05 2021

@author: dmattox
"""

import pickle
import umap
# import collections
import random
import Levenshtein

import plotnine as p9
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt


# Load encodings
with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/glycan_common_16048_may22.pickle', 'rb') as inFH:
    enc = pickle.load(inFH)

# Load SugarBase entries
with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/commonGlys_encodings.pickle', 'rb') as inFH:
    commonGlys = pickle.load(inFH)

# Load embeddings
with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/rep_withid_new.pkl', 'rb') as inFH:
    emb = pickle.load(inFH)

# Load embeddings from all glycans (train + validation)
with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/allsbid.pkl', 'rb') as inFH:
    allEmb = pickle.load(inFH)



emb['SBID4'].shape # Embedding from start token, represents entirity of the glycan
enc['SBID4']

enc['SBID4'][1]
enc['SBID4'][3]




allEmbs = np.array(list(emb.values()))

random.seed(27)


reducer = umap.UMAP().fit(allEmbs)

# with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/valGlyEmbUMAP_reducer.pkl', 'wb') as outFH:
#     pickle.dump(reducer, outFH, protocol=pickle.HIGHEST_PROTOCOL)

with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/valGlyEmbUMAP_reducer.pkl', 'rb') as inFH:
    reducer = pickle.load(inFH)

uu = reducer.embedding_

linkDict = {-1 : 'Unknown',
            0 : 'Free',
            1 : 'N-linked',
            2 : 'O-linked'}
linkTag = np.array([True if enc[k][1] != -1 else False for k in emb.keys()])
np.sum(linkTag)

immDict = {-1 : 'Unknown',
            0 : 'Non-immunogenic',
            1 : 'Immunogenic'}
immTag = np.array([True if enc[k][3] != -1 else False for k in emb.keys()])
np.sum(immTag)

kingDict = {-1 : 'Unknown',
            0 : 'Bacteria',
            1 : 'Metazoa',
            2 : 'Viridiplantae',
            3 : 'Fungi',
            4 : 'Virus',
            5 : 'Archaea'}
kingTag = np.array([True if -1 not in enc[k][2][1] else False for k in emb.keys()])
sum(kingTag)

uup = pd.DataFrame(uu, columns = ['UMAP1', 'UMAP2'])
uup['Linkage'] = [linkDict[enc[k][1]] for k in emb.keys()]
uup['Immunogenicity'] = [immDict[enc[k][3]] for k in emb.keys()]
uup['Kingdom'] = [kingDict[max(enc[k][2][1])] for k in emb.keys()] # If multiple labels, take the most rare one

plt = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2') + p9.geom_point(alpha = 0.2, color = 'grey')
plt

# plt + \
#     p9.scales.xlim(-6, 17.5) + p9.scales.ylim(-7, 10)


p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', color = 'Immunogenicity', alpha = 'Immunogenicity') + p9.geom_point(alpha = 0.2, color = 'grey') + \
    p9.geom_point() + \
        p9.scale_color_manual(values=['blue', 'orange', 'grey']) + p9.scale_alpha_manual(values = [0.6, 0.6, 0])

p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', color = 'Linkage', alpha = 'Linkage') + p9.geom_point(alpha = 0.2, color = 'grey') + \
    p9.geom_point() + \
        p9.scale_color_manual(values=['#cc9d1b', '#1bccc6', '#cc1b77', 'grey']) + p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0])

p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', color = 'Kingdom', alpha = 'Kingdom') + p9.geom_point(alpha = 0.2, color = 'grey') + \
    p9.geom_point() + \
        p9.scale_color_manual(["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "white", "#5153DE", "#9400D3"]) + \
            p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0.6, 0, 0.6, 0.6]) + \
                p9.theme_bw()

plt = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Kingdom', alpha = 'Kingdom') + p9.geom_point(alpha = 0.2, color = 'grey') + \
    p9.geom_point(size = 2) + \
        p9.scale_fill_manual(["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "white", "#5153DE", "#9400D3"]) + \
            p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0.6, 0, 0.6, 0.6]) + \
                p9.theme_minimal()

# plt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/glycan_embedding_kingdoms.pdf', dpi = 600)


ptcex = 3

tmpPltImm = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Kingdom', alpha = 'Kingdom') + p9.geom_point(alpha = 0.2, color = '#FFFFFF00', size = ptcex) + \
    p9.geom_point(size = ptcex, color = '#FFFFFF00') + \
        p9.scale_fill_manual(["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "#BEBEBE1A", "#5153DE", "#9400D3"]) + \
            p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0.6, 0, 0.6, 0.6]) + \
                p9.geom_point(uup[immTag], p9.aes(color = 'Immunogenicity'), fill = '#FFFFFF00', size = ptcex, alpha = 0.6) + \
                    p9.scale_color_manual(values = ['#a84c3e','#3e6ca8'])
tmpPltImm

plt_imm = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', fill = 'Kingdom', alpha = 'Kingdom') + p9.geom_point(alpha = 0.2, color = '#FFFFFF00', size = ptcex) + \
    p9.geom_point(size = ptcex, color = '#FFFFFF00') + \
        p9.scale_fill_manual(["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "#BEBEBE1A", "#5153DE", "#9400D3"]) + \
            p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0.6, 0, 0.6, 0.6]) + \
                p9.geom_point(uup[immTag], p9.aes(color = 'Immunogenicity'), fill = '#FFFFFF00', size = ptcex, alpha = 0.6) + \
                    p9.scale_color_manual(values = ['#a84c3e','#3e6ca8']) + \
                        p9.theme_minimal()

# plt_imm.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/glycan_embedding_king_imm.pdf', dpi = 600)

gDict = {g.sbID : g for g in commonGlys}

commonGlys = [g for g in commonGlys if g.sbID in emb.keys()]
gID2iupac = {g.sbID : g.iupac for g in commonGlys}


###
# Root Glc
'SBID4613' in gID2iupac.keys()
pnt = reducer.transform(emb['SBID4613'].reshape(1, -1))
tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)

'SBID' + '12054' in gID2iupac.keys()
pnt = reducer.transform(emb['SBID' + '12054'].reshape(1, -1))
tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)


'SBID' + '13622' in gID2iupac.keys() # 3Fuc, b4 Glc of 4613

#imm
'SBID' + '3980' in gID2iupac.keys() # 3Fuc, aGlc of 4613 & 13622
'SBID' + '1965' in gID2iupac.keys() # 2Fuc, aGlc of 12054


pnt = reducer.transform(allEmb['SBID' + '3980'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
tmpPlt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/tmp_glycan_embedding_imm.pdf', dpi = 600)

pnt = reducer.transform(allEmb['SBID' + '1965'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
tmpPlt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/tmp_glycan_embedding_imm.pdf', dpi = 600)

###
# Root GalNAc

'SBID' + '12833' in gID2iupac.keys() 
'SBID' + '5978' in gID2iupac.keys()

pnt1 = reducer.transform(allEmb['SBID' + '12833'].reshape(1, -1))[0]
pnt2 = reducer.transform(allEmb['SBID' + '5978'].reshape(1, -1))[0]

pnt1
tmpPltImm + p9.geom_point(p9.aes(x=pnt1[0], y=pnt1[1]), colour="magenta", size = 3)

pnt2
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt2[0], y=pnt2[1]), colour="magenta", size = 3)
tmpPlt
tmpPlt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/tmp_glycan_embedding_imm.pdf', dpi = 600)

gDict['SBID' + '2468'].print()
'SBID' + '2468' in gID2iupac.keys() # non-imm
pnt1 = reducer.transform(allEmb['SBID' + '2468'].reshape(1, -1))[0]
print(pnt1)
tmpPltImm + p9.geom_point(p9.aes(x=pnt1[0], y=pnt1[1]), colour="magenta", size = 3)




tmpPltImm + \
    p9.scales.ylim(7.5, 11) + p9.scales.xlim(0, 2)


immOClusTag = (uup['UMAP1'] > .75) & (uup['UMAP2'] > 8) & (uup['UMAP1'] < 1.25) & (uup['UMAP2'] < 9)
uup[immOClusTag]
[(g, gDict[g].immunogenic) for g,b in zip(emb.keys(), immOClusTag) if b]

gDict['SBID1573'].print()

tmpPltImm + \
    p9.scales.ylim(7.5, 11) + p9.scales.xlim(0, 2)

tmpPltImm + \
    p9.scales.ylim(8, 9) + p9.scales.xlim(.75, 1.25)



[g.print() for g in commonGlys if g.immunogenic == 'Yes' and g.link == 'N']

[g.print() for g in commonGlys if g.immunogenic == 'Yes' and g.link == 'O']

'SBID' + '12644' in gID2iupac.keys() # imm
pnt = reducer.transform(allEmb['SBID' + '12644'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt



'SBID' + '3897' in gID2iupac.keys() # imm
pnt = reducer.transform(allEmb['SBID' + '3897'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
tmpPlt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/tmp_glycan_embedding_imm.pdf', dpi = 600)

gDict['SBID' + '3897'].print()
gDict['SBID' + '3897'].taxonomy[1]

'SBID' + '1803' in gID2iupac.keys() # imm
pnt = reducer.transform(allEmb['SBID' + '1803'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
gDict['SBID' + '1803'].print()
gDict['SBID' + '1803'].taxonomy[1]



'SBID' + '15550' in gID2iupac.keys() # non-imm
pnt = reducer.transform(allEmb['SBID' + '15550'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
gDict['SBID' + '15550'].print()
gDict['SBID' + '15550'].taxonomy[1]


'SBID' + '13039' in gID2iupac.keys() # non-imm
pnt = reducer.transform(allEmb['SBID' + '13039'].reshape(1, -1))
tmpPlt = tmpPltImm + p9.geom_point(p9.aes(x=pnt[0][0], y=pnt[0][1]), colour="magenta", size = 3)
tmpPlt
gDict['SBID' + '13039'].print()
gDict['SBID' + '13039'].taxonomy[1]


# pairs = []
# uDists = []
# iSim = []

# embIDs = list(emb.keys())

# for i,gID1 in enumerate(embIDs):
#     if i % 50 == 0:
#         print(i)
#     for j in range(i, len(embIDs)):
#         if i != j:
#             gID2 = embIDs[j]
#             uDists.append(np.sqrt((uup.iloc[i]['UMAP1'] - uup.iloc[j]['UMAP1'])**2 + \
#                                   (uup.iloc[i]['UMAP2'] - uup.iloc[j]['UMAP2'])**2))
#             iSim.append(Levenshtein.ratio(gID2iupac[gID1], gID2iupac[gID2]))
#             pairs.append((gID1, gID2))

# plt.scatter(x = iSim, y = uDists, alpha = 0.1)

# gDict = {g.sbID : g for g in commonGlys}

# tmp = [p for i,u,p in zip(iSim, uDists, pairs) if i > 0.85 and 10 < u < 20]
# for p in tmp:
#     gDict[p[0]].print()
#     print()
#     gDict[p[1]].print()
#     print('---------------')

# # Try tsne

# import sys; sys.path.append('/Users/dmattox/FIt-SNE')
# from fast_tsne import fast_tsne
# from sklearn.decomposition import PCA

# pca = PCA()
# pca.fit(np.transpose(allEmbs))

# pca.components_.shape

# pca.components_[0,:].shape
# plt.scatter(x = pca.components_[0,:], y = pca.components_[1,:])

# X = np.transpose(pca.components_[0:50,:])

# %time tsne_default = fast_tsne(X, seed=42)

# tsp = pd.DataFrame(tsne_default, columns = ['tSNE1', 'tSNE2'])
# tsp['Linkage'] = [linkDict[enc[k][1]] for k in emb.keys()]
# tsp['Immunogenicity'] = [immDict[enc[k][3]] for k in emb.keys()]
# tsp['Kingdom'] = [kingDict[max(enc[k][2][1])] for k in emb.keys()] # If multiple labels, take the most rare one

# p9.ggplot(tsp) + p9.aes(x = 'tSNE1', y = 'tSNE2', color = 'Immunogenicity', alpha = 'Immunogenicity') + p9.geom_point(alpha = 0.2, color = 'grey') + \
#     p9.geom_point() + \
#         p9.scale_color_manual(values=['blue', 'orange', 'grey']) + p9.scale_alpha_manual(values = [0.6, 0.6, 0])
    
# p9.ggplot(tsp) + p9.aes(x = 'tSNE1', y = 'tSNE2', color = 'Linkage', alpha = 'Linkage') + p9.geom_point(alpha = 0.2, color = 'grey') + \
#     p9.geom_point() + \
#         p9.scale_color_manual(values=['#cc9d1b', '#1bccc6', '#cc1b77', 'grey']) + p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0])
    

# p9.ggplot(tsp) + p9.aes(x = 'tSNE1', y = 'tSNE2', color = 'Kingdom', alpha = 'Kingdom') + p9.geom_point(alpha = 0.2, color = 'grey') + \
#     p9.geom_point() + \
#         p9.scale_color_manual(["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "white", "#5153DE", "#9400D3"]) + \
#             p9.scale_alpha_manual(values = [0.6, 0.6, 0.6, 0.6, 0, 0.6, 0.6]) + \
#                 p9.theme_bw()




















