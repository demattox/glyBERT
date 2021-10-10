#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:01:25 2021

@author: dmattox
"""

import pickle, umap, collections, random
import matplotlib.pyplot

import plotnine as p9

import numpy as np
import pandas as pd

# import matplotlib.pyplot as plt

# from plotnine import ggplot, aes, geom_line


# emb = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/emb_label2.pkl', 'rb'))
# emb = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/emb_label_withid.pkl', 'rb'))
emb = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/emb_label_withid_new.pkl', 'rb'))

# Load SugarBase entries
with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/commonGlys_encodings.pickle', 'rb') as inFH:
    commonGlys = pickle.load(inFH)
gDict = {g.sbID : g for g in commonGlys}

monoDict = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/commonMonoDict.pickle', 'rb'))
revMonoDict = {v:k for k,v in monoDict.items()}

len(emb)
len(emb[0])
emb[0]['SBID5000'].shape
emb[0]['SBID5000'][0].shape
len(emb[1])
emb[1]['SBID5000'].shape

# emb[1] = np.array([revMonoDict[int(m) - 2] for m in emb[1]])

allMonoLabels = [ revMonoDict[int(num) - 2] for v in emb[1].values() for num in v ]
allEmbs = np.array([ arr for v in emb[0].values() for arr in v])
allEmbs.shape


monoCnts = collections.Counter(allMonoLabels)
monoCnts.most_common()

topMono = np.array([k for k,v in monoCnts.items() if v > 100])

topTag = np.array([True if e in topMono else False for e in allMonoLabels])


random.seed(18)

reducer = umap.UMAP().fit(allEmbs)
# uu = reducer.fit_transform(allEmbs)


# with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/valMonoEmbUMAP_reducer.pkl', 'wb') as outFH:
#     pickle.dump(reducer, outFH, protocol=pickle.HIGHEST_PROTOCOL)

with open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/valMonoEmbUMAP_reducer.pkl', 'rb') as inFH:
    reducer = pickle.load(inFH)


uu = reducer.embedding_


# plt.scatter(x = uu[:,0], y = uu[:,1], alpha = 0.0055, c = 'grey')

# plt.figure(figsize=(15,15))
# plt.scatter(x = uu[:,0], y = uu[:,1], alpha = 0.0055, c = 'grey')
# plt.scatter(x = uu[topTag,0], y = uu[topTag,1], alpha = 0.0055, c = 'red')



uup = pd.DataFrame(uu, columns = ['UMAP1', 'UMAP2'])



uup['Monosacc'] = [l if l in topMono else 'other' for l in allMonoLabels]

uup['All M']

plt = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', color = 'Monosacc') + p9.geom_point(alpha = 0.5)

plt = plt +  p9.scale_color_manual(values=[#"#ff0000", # Fuc/red
                                   "#00A651", #Rha/green
                                   "#ff0000", # Fuc/red
                                   "#ffd900", # Gal/yellow
                                   "#ffd900", 
                                   "#ffd900",
                                   "#ffd900",
                                   "#0080ff", # Glc/blue
                                   "#0080ff",
                                   "#0080ff",
                                   "#0080ff",
                                   "#ffd900", # kdo, yellow
                                   "#00A651", # LDManHep, green
                                   "#00A651", # Man
                                   "#9e1fff", # NeuAc, purple
                                   "#96f2f7", # NeuGc, light blue
                                   "#00A651", #Rha/green
                                   "#ff5900", # Xyl, orange
                                   "grey"]) # other

plt

plt_leg = p9.ggplot(uup) + p9.aes(x = 'UMAP1', y = 'UMAP2', color = 'Monosacc') + p9.geom_point(alpha = 1) +  p9.scale_color_manual(values=[#"#ff0000", # Fuc/red
                                   "#00A651", #Rha/green
                                   "#ff0000", # Fuc/red
                                   "#ffd900", # Gal/yellow
                                   "#ffd900", 
                                   "#ffd900",
                                   "#ffd900",
                                   "#0080ff", # Glc/blue
                                   "#0080ff",
                                   "#0080ff",
                                   "#0080ff",
                                   "#ffd900", # kdo, yellow
                                   "#00A651", # LDManHep, green
                                   "#00A651", # Man
                                   "#9e1fff", # NeuAc, purple
                                   "#96f2f7", # NeuGc, light blue
                                   "#00A651", #Rha/green
                                   "#ff5900", # Xyl, orange
                                   "grey"])  + p9.theme_minimal() # other 

plt_leg

# plt_leg.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/monosacc_embedding-legend.pdf', dpi = 600)


labCoords = []
for sugar in topMono:
    labCoords.append(np.mean(np.array([c for l,c in zip(allMonoLabels, uu) if l == sugar]), axis = 0))
labCoords = np.array(labCoords)

mean_uup = pd.DataFrame()

mean_uup['Monosacc'] = topMono
mean_uup['UMAP1'] = labCoords[:,0]
mean_uup['UMAP2'] = labCoords[:,1]

plt = plt + p9.geom_text(data = mean_uup, label = mean_uup['Monosacc'], color = 'black') + p9.theme_minimal()

plt

# plt.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/monosacc_embedding.pdf', dpi = 600)1


enc = pickle.load(open('/Users/dmattox/Documents/qbs/cbk/sugarTrees/pickles/glycan_common_16048_may22.pickle', 'rb'))

enc_val = {k:v for k,v in enc.items() if k in emb[0]} # All encodings for glycans in validatio

len(enc_val)
len(emb[0])

coreFuc_glys = []
for gly in enc_val.keys():
    if monoDict['Fuc'] in enc_val[gly][0][:,0]:
        fTag = np.logical_and(enc_val[gly][0][:,0] == monoDict['Fuc'], enc_val[gly][0][:,5] == 2) # true for row in array with Fucose at depth two (if present)
        if sum(fTag) == 1: # if a fuc residue present at depth 2
            if enc_val[gly][0][fTag][0,4] == 6 and enc_val[gly][0][enc_val[gly][0][:,5] == 1][0,0] == monoDict['GlcNAc']:# alpha 1-6 link to GlcNAc (core fucosylation)
                coreFuc_glys.append(gly)

enc_val['SBID19184']
emb[1]['SBID19184'] == (monoDict['Fuc'] + 2)

coreFuc_embs = []
for gly in coreFuc_glys:
    coreFuc_embs.append(emb[0][gly][emb[1][gly] == (monoDict['Fuc'] + 2)][-1]) # Get all embeddings for fucose monosaccharides from glycans with a terminal fucose, then take the last one in the array (where core fucsylations end up in the encoding)

len(coreFuc_embs)
len(coreFuc_glys)

coreFuc_embs = np.array(coreFuc_embs)

uu_cF = reducer.transform(coreFuc_embs)

np.mean(uu_cF, axis = 0)
mean_uucF = pd.DataFrame()
mean_uucF['UMAP1'] = [np.mean(uu_cF, axis = 0)[0]]
mean_uucF['UMAP2'] = [np.mean(uu_cF, axis = 0)[1]]
mean_uucF['Monosacc'] = 'Fuc'

uu_cF = pd.DataFrame(uu_cF, columns = ['UMAP1', 'UMAP2'])
uu_cF['Monosacc'] = 'Fuc'

plt_cF = plt + p9.geom_point(data = uu_cF, color = 'black', shape = '<', alpha = 0.01) + \
    p9.geom_label(data = mean_uucF, label = 'Core\nFuc', color = 'black', alpha = 0.3)


plt_cF

# plt_cF.save('/Users/dmattox/Documents/qbs/cbk/sugarTrees/plots/monosacc_embedding_coreFuc.pdf', dpi = 600)

np.max(uup['UMAP2'][uup['Monosacc'] == 'Fuc'])


# Nonulsonic acid cluster

neuClusTag = (uup['UMAP1'] > -5) & (uup['UMAP2'] > 18.5) & (uup['UMAP1'] < 1)
collections.Counter([m for m,b in zip(allMonoLabels, neuClusTag) if b]).most_common()

# Root Glc

rootGlc_glys = []
for gly in enc_val.keys():
    gTag = np.logical_and(enc_val[gly][0][:,0] == monoDict['Glc'], enc_val[gly][0][:,5] == 1) # true for row in array with Glc at depth 1 (if present)
    if sum(gTag) == 1: # if a glc residue present at depth 1
        rootGlc_glys.append(gly)

len(rootGlc_glys)

gly = rootGlc_glys[0]
all([enc_val[gly][0][1,0] == 4 for gly in rootGlc_glys])

rootGlc_embs = []
for gly in rootGlc_glys:
    rootGlc_embs.append(emb[0][gly][emb[1][gly] == (monoDict['Glc'] + 2)][0]) # Get all embeddings for the first Glc monosaccharides from glycans with a root glucose

len(rootGlc_glys)
len(rootGlc_embs)

rootGlc_embs = np.array(rootGlc_embs)

uu_rG = reducer.transform(rootGlc_embs)


np.mean(uu_rG, axis = 0)
mean_uurG = pd.DataFrame()
mean_uurG['UMAP1'] = [np.mean(uu_rG, axis = 0)[0]]
mean_uurG['UMAP2'] = [np.mean(uu_rG, axis = 0)[1]]
mean_uurG['Monosacc'] = 'Glc'

uu_rG = pd.DataFrame(uu_rG, columns = ['UMAP1', 'UMAP2'])
uu_rG['Monosacc'] = 'Glc'

plt_rG = plt + p9.geom_point(data = uu_rG, color = 'black', shape = '<', alpha = 0.05) + \
    p9.geom_label(data = mean_uurG, label = 'Root\nGlc', color = 'black', alpha = 0.3)


plt_rG

collections.Counter([gDict[gly].link for gly in rootGlc_glys])

plt_rG + p9.scales.xlim(-7, -4) + p9.scales.ylim(3.5, 7)

rootGlcClusTag = (uup['UMAP1'] > -7) & (uup['UMAP2'] > 4) & (uup['UMAP1'] < -4) & (uup['UMAP2'] < 5.75)
collections.Counter([m for m,b in zip(allMonoLabels, rootGlcClusTag) if b]).most_common()




# Root GalNAc
rootGalNAc_glys = []
for gly in enc_val.keys():
    gTag = np.logical_and(enc_val[gly][0][:,0] == monoDict['GalNAc'], enc_val[gly][0][:,5] == 1) # true for row in array with Glc at depth 1 (if present)
    if sum(gTag) == 1: # if a glc residue present at depth 1
        rootGalNAc_glys.append(gly)

len(rootGalNAc_glys)

gly = rootGalNAc_glys[0]
all([enc_val[gly][0][1,0] == monoDict['GalNAc'] for gly in rootGalNAc_glys])

rootGalNAc_embs = []
for gly in rootGalNAc_glys:
    rootGalNAc_embs.append(emb[0][gly][emb[1][gly] == (monoDict['GalNAc'] + 2)][0]) # Get all embeddings for the first Glc monosaccharides from glycans with a root glucose

len(rootGalNAc_glys)
len(rootGalNAc_embs)

rootGalNAc_embs = np.array(rootGalNAc_embs)

uu_rG = reducer.transform(rootGalNAc_embs)

np.mean(uu_rG, axis = 0)
mean_uurG = pd.DataFrame()
mean_uurG['UMAP1'] = [np.mean(uu_rG, axis = 0)[0]]
mean_uurG['UMAP2'] = [np.mean(uu_rG, axis = 0)[1]]
mean_uurG['Monosacc'] = 'GalNAc'

uu_rG = pd.DataFrame(uu_rG, columns = ['UMAP1', 'UMAP2'])
uu_rG['Monosacc'] = 'GalNAc'

plt_rG = plt + p9.geom_point(data = uu_rG, color = 'black', shape = '<', alpha = 0.01) + \
    p9.geom_label(data = mean_uurG, label = 'Root\nGalNAc', color = 'black', alpha = 0.3)


lGlys = []
rGlys = []
for i,gly in enumerate(rootGalNAc_glys):
    if uu_rG['UMAP1'][i] < 0:
        lGlys.append(gly)
    else:
        rGlys.append(gly)

len(lGlys)
len(rGlys)


collections.Counter([gDict[g].link for g in lGlys]).most_common()
collections.Counter([gDict[g].immunogenic for g in lGlys]).most_common()
# [gDict[g].iupac for g in lGlys]
# [gDict[g].taxonomy[0] for g in lGlys if gDict[g].taxonomy]
# [gDict[g].species for g in lGlys if gDict[g].taxonomy]
collections.Counter([gDict[g].taxonomy[0][0] for g in lGlys if gDict[g].taxonomy]).most_common()
collections.Counter(['-'.join(list(gDict[g].tree[(1,0)].Clinks.keys())) for g in lGlys]).most_common()
# [[v.base for n,v in gDict[g].tree.items() if v.depth == 2] for g in lGlys]
collections.Counter([gDict[g].tree[(1,0)].anomeric for g in lGlys]).most_common()
collections.Counter([max(gDict[g].tree[(1,0)].subway) for g in lGlys]).most_common()
collections.Counter([len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in lGlys]).most_common()
matplotlib.pyplot.hist([len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in lGlys])

yl = [max(gDict[g].tree[(1,0)].subway) for g in lGlys]
xl = [len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in lGlys]
cl = ['blue' for g in lGlys]




collections.Counter([gDict[g].link for g in rGlys]).most_common()
collections.Counter([gDict[g].immunogenic for g in rGlys]).most_common()
# [gDict[g].iupac for g in rGlys]
# [gDict[g].taxonomy[0] for g in rGlys if gDict[g].taxonomy]
collections.Counter([gDict[g].taxonomy[0][0] for g in rGlys if gDict[g].taxonomy]).most_common()
collections.Counter(['-'.join(list(gDict[g].tree[(1,0)].Clinks.keys())) for g in rGlys]).most_common()
# [[v.base for n,v in gDict[g].tree.items() if v.depth == 2] for g in rGlys]
collections.Counter([gDict[g].tree[(1,0)].anomeric for g in rGlys]).most_common()
collections.Counter([max(gDict[g].tree[(1,0)].subway) for g in rGlys]).most_common()
collections.Counter([len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in rGlys]).most_common()
matplotlib.pyplot.hist([len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in rGlys])

yr = [max(gDict[g].tree[(1,0)].subway) for g in rGlys]
xr = [len([n for n,v in gDict[g].tree.items() if v.base not in ['START', 'END']]) for g in rGlys]
cr = ['red' for g in lGlys]

matplotlib.pyplot.hist2d(x = xl + xr, y = yl + yr)

matplotlib.pyplot.hist2d(x = xl, y = yl, bins = 5)
matplotlib.pyplot.ylabel('Number of branches')
matplotlib.pyplot.xlabel('Number of sacc. residues')

matplotlib.pyplot.hist2d(x = xr, y = yr, bins = 5)
matplotlib.pyplot.ylabel('Number of branches')
matplotlib.pyplot.xlabel('Number of sacc. residues')

###

tmp = [g for g in gDict.values() if g.immunogenic == 'Yes' and g.link == 'O']
len(tmp)
for g in tmp:
    g.print()

rootGlys = [g.tree[(1,0)].base for g in gDict.values()]
len(rootGlys)
collections.Counter(rootGlys)


# Root GlcNAc
rootGlcNAc_glys = []
for gly in enc_val.keys():
    gTag = np.logical_and(enc_val[gly][0][:,0] == monoDict['GlcNAc'], enc_val[gly][0][:,5] == 1) # true for row in array with Glc at depth 1 (if present)
    if sum(gTag) == 1: # if a glc residue present at depth 1
        rootGlcNAc_glys.append(gly)

len(rootGlcNAc_glys)

gly = rootGlcNAc_glys[0]
all([enc_val[gly][0][1,0] == monoDict['GlcNAc'] for gly in rootGlcNAc_glys])

rootGlcNAc_embs = []
for gly in rootGlcNAc_glys:
    rootGlcNAc_embs.append(emb[0][gly][emb[1][gly] == (monoDict['GlcNAc'] + 2)][0]) # Get all embeddings for the first Glc monosaccharides from glycans with a root glucose

len(rootGlcNAc_glys)
len(rootGlcNAc_embs)

rootGlcNAc_embs = np.array(rootGlcNAc_embs)

uu_rG = reducer.transform(rootGlcNAc_embs)

np.mean(uu_rG, axis = 0)
mean_uurG = pd.DataFrame()
mean_uurG['UMAP1'] = [np.mean(uu_rG, axis = 0)[0]]
mean_uurG['UMAP2'] = [np.mean(uu_rG, axis = 0)[1]]
mean_uurG['Monosacc'] = 'GlcNAc'

uu_rG = pd.DataFrame(uu_rG, columns = ['UMAP1', 'UMAP2'])
uu_rG['Monosacc'] = 'GlcNAc'

plt_rG = plt + p9.geom_point(data = uu_rG, color = 'black', shape = '<', alpha = 0.1) + \
    p9.geom_label(data = mean_uurG, label = 'Root\nGlcNAc', color = 'black', alpha = 0.3)

plt_rG

rGlys = []
lGlys = []
for i,gly in enumerate(rootGlcNAc_glys):
    if uu_rG['UMAP1'][i] < 0:
        lGlys.append(gly)
    else:
        rGlys.append(gly)

len(lGlys)
len(rGlys)



collections.Counter([gDict[g].link for g in lGlys]).most_common()
collections.Counter([gDict[g].immunogenic for g in lGlys]).most_common()
# [gDict[g].iupac for g in lGlys]
# [gDict[g].taxonomy[0] for g in lGlys if gDict[g].taxonomy]
# [gDict[g].species for g in lGlys if gDict[g].taxonomy]
collections.Counter([gDict[g].taxonomy[0][0] for g in lGlys if gDict[g].taxonomy]).most_common()
collections.Counter(['-'.join(list(gDict[g].tree[(1,0)].Clinks.keys())) for g in lGlys]).most_common()
# [[v.base for n,v in gDict[g].tree.items() if v.depth == 2] for g in lGlys]
collections.Counter([gDict[g].tree[(1,0)].anomeric for g in lGlys]).most_common()
collections.Counter([max(gDict[g].tree[(1,0)].subway) for g in lGlys]).most_common()
collections.Counter([len(gDict[g].tree) for g in lGlys]).most_common()
matplotlib.pyplot.hist([len(gDict[g].tree) for g in lGlys])


collections.Counter([gDict[g].link for g in rGlys]).most_common()
collections.Counter([gDict[g].immunogenic for g in rGlys]).most_common()
# [gDict[g].iupac for g in rGlys]
# [gDict[g].taxonomy[0] for g in rGlys if gDict[g].taxonomy]
collections.Counter([gDict[g].taxonomy[0][0] for g in rGlys if gDict[g].taxonomy]).most_common()
collections.Counter(['-'.join(list(gDict[g].tree[(1,0)].Clinks.keys())) for g in rGlys]).most_common()
# [[v.base for n,v in gDict[g].tree.items() if v.depth == 2] for g in rGlys]
collections.Counter([gDict[g].tree[(1,0)].anomeric for g in rGlys]).most_common()
collections.Counter([max(gDict[g].tree[(1,0)].subway) for g in rGlys]).most_common()
collections.Counter([len(gDict[g].tree) for g in rGlys]).most_common()
matplotlib.pyplot.hist([len(gDict[g].tree) for g in rGlys])

