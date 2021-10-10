#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 08:23:57 2021

@author: dmattox
"""

import pickle, collections
import numpy as np
import matplotlib.pyplot as plt

import cbk.sugarTrees.bin.sugarTrees as SugarTrees

homeDir = '/Users/dmattox/cbk/sugarTrees/'

sbFile = homeDir+'data/SugarBase_v2_2020_11_19.csv'

sbEntries = []

with open(sbFile, 'r') as inFH: # Read in sugar base csv file
    for line in inFH:
        line = line.strip()
        line = line.split('","')
        if line[0] != '"SugarBase_ID':
            line = [l.replace('"','') for l in line] # Quoted csv, strip front quote not caught with split
            sbEntries.append(SugarTrees.SugarBase(*line)) # append database entry to list as SugarBase object

okayEntries = []
failed = []
failedIUPACs = []

for i, sugar in enumerate(sbEntries):
    try:
        sugar.buildTree()
        okayEntries.append(sugar)
    except:
        #print('FAILED', sugar.sbID)
        # if sugar.iupac[-1] == ')':
        #     print(sugar.sbID, sugar.iupac)
        failed.append(sugar.sbID)
        failedIUPACs.append(sugar.iupac)

print(len(okayEntries), 'correctly formatted glycans') # 18154
len(failed) # 1145

len(okayEntries) + len(failed) #19299
len(sbEntries)

len(failed)/len(sbEntries)

for sugar in okayEntries:
    sugar.drawSbwyMap()

#######################
# Duplicated glycans?
#######################

uniGlyIds = collections.defaultdict(list)
for gly in okayEntries:
    gID = list(gly.tree.keys())
    gID.extend([tuple(v.Clinks.keys()) for v in gly.tree.values()]) # Add tuples containing linkage info to the tree structure key
    gID.append(''.join([node.anomeric for node in gly.tree.values()])) # Add string for unique anomeric conformations to the tree structure key
    gID.append('-'.join([node.base for node in gly.tree.values()])) # Add concatenated monosaccaride names
    uniGlyIds[tuple(gID)].append(gly)

print(len(uniGlyIds), 'unique glycans') # 16,835 unique glycans
len(okayEntries) - len(uniGlyIds) # 1,319 replicated glycans

duplicatedGlys = {k:v for k,v in uniGlyIds.items() if len(v) > 1}
len(duplicatedGlys)

collections.Counter([len(v) for v in uniGlyIds.values()]).most_common()

manyDups = [v for v in uniGlyIds.values() if len(v) > 3]

# for g in manyDups[2]:
#     # print(g.iupac)
#     g.print()
#     g.treePrint()

# Check examples
ex1 = [g for g in okayEntries if g.sbID == 'SBID8049'][0]
ex2 = [g for g in okayEntries if g.sbID == 'SBID12263'][0]
ex1.iupac
ex2.iupac

ex1 = [g for g in okayEntries if g.sbID == 'SBID8088'][0]
ex2 = [g for g in okayEntries if g.sbID == 'SBID13110'][0]
ex1.iupac
ex2.iupac


for repList in duplicatedGlys.values():
    glySizes = [len(gly.tree) for gly in repList]
    if len(set(glySizes)) != 1:
        print(repList)


# How are duplicated glycans labelled?
dupImm = []
dupLink = []
dupSpec = []
for dupGlys in duplicatedGlys.values():
    dupImm.append(len(set([g.immunogenic for g in dupGlys if g.immunogenic != 'Unknown'])))
    dupLink.append(len(set([g.link for g in dupGlys if g.link != 'None'])))
    dupSpec.append(len(set([tuple(g.species) for g in dupGlys if g.species != ['']])))
    # if len(set([tuple(g.species) for g in dupGlys if g.species != ['']])) > 1:
    #     print(set([tuple(g.species) for g in dupGlys if g.species != ['']]))

collections.Counter(dupImm)
collections.Counter(dupLink)
collections.Counter(dupSpec)

# If one replicate is labelled and other isn't, use known label
# If replicates have conflicting labels, set label to unknwon
# If replicates have different species of origin, take the set of reported species

# Merge/reconcile labels for the lowest numbered SBID, drop replicates
delLst = [] # Initialize list to record duplicated/replicated glycans to remove

for dupGlys in duplicatedGlys.values():
    allLink = set(g.link for g in dupGlys)
    allImm = set(g.immunogenic for g in dupGlys)
    allSpec = set(s for g in dupGlys for s in g.species)
    
    recGly = dupGlys[np.argmin(g.id for g in dupGlys)] # Update labels for the glycan with the lowest SBID #, modifies corresponding glycan in okayEntries list
    delLst.extend([g for g in dupGlys if g != recGly]) # Mark the other ones to remove from the okayEntries list
    
    # Immunogenicity
    if ('No' in allImm) or ('Yes' in allImm): # If any glycans are labelled
        if ('No' in allImm) and ('Yes' in allImm): # if labels conflict
            recGly.immunogenic = 'Unknown'
        else:
            recGly.immunogenic = [l for l in allImm if l != 'Unknown'][0]
    else:
        recGly.immunogenic = 'Unknown'
        
    # Linkage
    infoLinks = [l for l in allLink if l != 'None']
    if len(infoLinks) > 1:
        recGly.link = 'None'
    elif len(infoLinks) == 1:
        recGly.link = infoLinks[0]
    else:
        recGly.link = 'None'
    
    # Species
    if allSpec == {''}:
        recGly.species = list(allSpec)
    else:
        recGly.species = [s for s in allSpec if s != '']

    
    

# [g for g in okayEntries if g.sbID == 'SBID18813'][0].print()
# [g for g in okayEntries if g.sbID == 'SBID9500'][0].print()

# Delete redundant glycans
len(okayEntries)
delLst = [i for i in range(len(okayEntries)) if okayEntries[i] in delLst] # Get glycan indices in okayEntries for replicated glys to delete
delLst.sort(reverse=True)
for i in delLst: # Loop through in reverse order to preserve indexing
    del okayEntries[i]
len(okayEntries)

###################
# Get IUPAC from encoding test

# sugar.print()
# sugar.treePrint()
# for node in sugar.tree.values():
#     print(node)
#     print(node.Clinks)

# sugar.buildEncoding(maxB = maxB, maxC = maxC, monoDict = monoDict)
# sugar.encoding

###################

# Read in genus to superkingdom mappings
g2rank_file = homeDir + 'pickles/genera2ranks.pickle'
with open(g2rank_file, 'rb') as inFH:
    taxa_relations = pickle.load(inFH) # list of dictionaries mapping genus to 0: superkingdom, 1: kingdom, 2: phylum, 3: class, 4: order, 5: family


# Get glycan stats
monosaccLst = set()
for sugar in okayEntries:
    monosaccLst.update([s.base for s in sugar.tree.values()])

monosaccLst = list(monosaccLst)
len(monosaccLst)

maxC = []
for sugar in okayEntries:
    try:
        maxC.append(int(max([max(s.Clinks.keys()) for s in sugar.tree.values()])))
        if int(max([max(s.Clinks.keys()) for s in sugar.tree.values()])) == 0:
            sugar.print()
            print()
    except:
        failed.append(sugar)
        okayEntries.remove(sugar)
        
cntC = collections.Counter(maxC)
cntC.most_common()
maxC = max(maxC)


maxB = [max(sugar.tree[(0,0)].subway) for sugar in okayEntries]

cntB = collections.Counter(maxB)
cntB.most_common()

sum(cntB.values())
len(okayEntries)

# Drop the outlier glycan with 24 branches
okayEntries = [sugar for sugar in okayEntries if max(sugar.tree[(0,0)].subway) != 24]
len(okayEntries)

maxB = [max(sugar.tree[(0,0)].subway) for sugar in okayEntries]

cntB = collections.Counter(maxB)
cntB.most_common()

maxB = max(maxB)


########
# Repeat monosaccLst with final set of 16,834 glycans
monosaccLst = set()
for sugar in okayEntries:
    monosaccLst.update([s.base for s in sugar.tree.values()])

monosaccLst = list(monosaccLst)
len(monosaccLst) # 944 unique monosaccharides

cntM = collections.Counter()
for sugar in okayEntries:
    shortLst = set([s.base for s in sugar.tree.values()])
    for m in shortLst:
        cntM[m] += 1
    
cntM.most_common()

Ms = [cnt[0] for cnt in cntM.most_common()]
cnts = [cnt[1] for cnt in cntM.most_common()]

cntCnts = collections.Counter(cnts)
sum([cntCnts[i] for i in [1,2,3]])

len(monosaccLst) - sum([cntCnts[i] for i in [1,2,3]])

if __name__ == 'main':
    plt.plot(list(range(len(cnts))), np.log10(cnts))
    plt.plot(list(range(len(cnts))), [0]*len(cnts), '--k')
    plt.vlines(415, 0, np.log10(3), colors = 'red')
    plt.show()

commonMs = Ms[:415]
commonGlys = []
for sugar in okayEntries:
    if all(s.base in commonMs  for s in sugar.tree.values()):
        commonGlys.append(sugar)
print(str(len(commonGlys)) + ' commonGlys') # If we limit to glycans containing the 413 most common monosaccharides [ + 2 start/stop tokens = 415 ] (appearing in > 3 glycans), 16048 (of 16834) glycans would remain, excludes 786 glycans
len(okayEntries) - len(commonGlys)



#############################

# monosacc dictionary (id --> index)
ordered_monosaccLst = commonMs[:] # List to hold ordered monosaccharides
# Add start and end tok//ens
ordered_monosaccLst.insert(0,ordered_monosaccLst.pop(ordered_monosaccLst.index('START')))
ordered_monosaccLst.insert(1,ordered_monosaccLst.pop(ordered_monosaccLst.index('END')))

monoDict = {}
for i, base in enumerate(ordered_monosaccLst):
    monoDict[base] = i



#anomeric specification dictionary
ano = []
for sugar in okayEntries:
    ano.extend(list(set([b.anomeric for b in sugar.tree.values()])))
    
set(ano)

collections.Counter(ano)
anoDict = {'u': 0, 'a': 1, 'b': 2}


# Check immunogenicity info
imm = []
for sugar in commonGlys:
    imm.append(sugar.immunogenic.strip())
# print(imm)
collections.Counter(imm)
collections.Counter(imm)['No'] + collections.Counter(imm)['Yes']

immDict = {'Unknown': -1, 'No': 0, 'Yes': 1}

# check links
lis = []
for sugar in commonGlys:
    lis.append(sugar.link.strip())
# print(lis)
collections.Counter(lis)


linkDict = {'None': -1, 'Free': 0, 'N': 1, 'O': 2}


#Species encoding
spec = set()
for sugar in okayEntries:
    spec.update([s.strip() for s in sugar.species])
# print(len(spec)) # 1538 unique species
spec = sorted(list(spec))

speciesDict = {i:s for i,s in enumerate(spec)}
speciesDictRev = {s:i for i,s in enumerate(spec)}

# Fill in sugar.taxonomy with full taxanomic information gathered from generas  (first 6 ranks)
for i,rank in enumerate(taxa_relations):
    for sugar in okayEntries:
        if sugar.species != ['']:
            taxa = set()
            for s in sugar.species:
                if s.split()[0] in rank.keys():
                    taxa.add(rank[s.split()[0]])
                else:
                    if s[-5:] == 'Virus' and i in [0,1]:  # if considering superkingdoms or kingdoms (1st/2nd rank in list of taxa relations), include viruses as separate taxa rank label
                        taxa.add('Virus')
                    # else: # no taxon is provided for this species at this rank
                        # print(s) 
            sugar.taxonomy.append(list(taxa)) # Add unique list of taxa for the given rank associated with the listed species
            
# # add genera and species to taxonomy list
# for sugar in okayEntries:
#     if sugar.species != ['']:
#         generaLst = set()
#         specLst = set()
#         for s in sugar.species:
            

# stats at each level of taxonomy across 16,048 glycans with more common monosaccharides
tax_stats = []
for i in range(6):
    tax_stats.append(collections.Counter()) # initialize list of counters, one for each level of taxonomy

for sugar in commonGlys:
    if sugar.species != ['']:
        for i in range(len(sugar.taxonomy)):
            for taxid in sugar.taxonomy[i]:
                tax_stats[i][taxid] += 1
                
orderedCounts = [list(zip(*t.most_common())) for t in tax_stats]


if __name__ == 'main':
    plt.bar(orderedCounts[0][0], orderedCounts[0][1])
    plt.bar(orderedCounts[1][0], orderedCounts[1][1])
    plt.bar(orderedCounts[2][0], orderedCounts[2][1])
    plt.bar(orderedCounts[3][0], orderedCounts[3][1])
    plt.bar(orderedCounts[4][0], orderedCounts[4][1])
    plt.bar(orderedCounts[5][0], orderedCounts[5][1])

# Stats at each level of taxonomy across all properly formatted glycans
tax_stats = []
for i in range(6):
    tax_stats.append(collections.Counter()) # initialize list of counters, one for each level of taxonomy

for sugar in okayEntries:
    if sugar.species != ['']:
        for i in range(len(sugar.taxonomy)):
            for taxid in sugar.taxonomy[i]:
                tax_stats[i][taxid] += 1
                
orderedCounts = [list(zip(*t.most_common())) for t in tax_stats]

if __name__ == 'main':
    plt.bar(orderedCounts[0][0], orderedCounts[0][1])
    plt.bar(orderedCounts[1][0], orderedCounts[1][1])
    plt.bar(orderedCounts[2][0], orderedCounts[2][1])
    plt.bar(orderedCounts[3][0], orderedCounts[3][1])
    plt.bar(orderedCounts[4][0], orderedCounts[4][1])
    plt.bar(orderedCounts[5][0], orderedCounts[5][1])



taxDicts = [] # holds dictionaries to map from taxonomic identifier to a number representation (sorted by frequency)
for i in range(6):
    tax = {'': -1}
    for j,taxid in enumerate(orderedCounts[i][0]):
        tax[taxid] = j
    taxDicts.append(tax)
    
# Check multilabel occurences in the common monosaccharide-containing glycans
multilabelStats = []
for i in range(6):
    multilabelStats.append(collections.Counter())
    
for sugar in commonGlys:
    if sugar.species != ['']:
        for i in range(len(sugar.taxonomy)):
            multilabelStats[i][len(sugar.taxonomy[i])] += 1

mLabelOrdered = [sorted(t.most_common()) for t in multilabelStats]
mLabelOrdered = [list(zip(*t)) for t in mLabelOrdered]

if __name__ == 'main':
    plt.bar(mLabelOrdered[0][0], mLabelOrdered[0][1])
    plt.title('superkingdoms')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()
    
    plt.bar(mLabelOrdered[1][0], mLabelOrdered[1][1])
    plt.title('kingdoms')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()
    
    plt.bar(mLabelOrdered[2][0], mLabelOrdered[2][1])
    plt.title('phyla')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()
    
    plt.bar(mLabelOrdered[3][0], mLabelOrdered[3][1])
    plt.title('classes')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()
    
    plt.bar(mLabelOrdered[4][0], mLabelOrdered[4][1])
    plt.title('orders')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()
    
    plt.bar(mLabelOrdered[5][0], mLabelOrdered[5][1])
    plt.title('families')
    plt.xlabel('number of labels/glycan')
    plt.ylabel('count')
    plt.show()