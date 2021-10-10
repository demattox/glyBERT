#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 08:15:32 2020

@author: dmattox
"""

from ete3 import NCBITaxa
import pickle, collections

ncbi = NCBITaxa()


with open('/Users/dmattox/cbk/sugarTrees/pickles/speciesDict.pickle', 'rb') as inFH:
    specDict = pickle.load(inFH)
    
allSpecs = list(specDict.values())
allSpecs.remove('')

genera = list(set([s.split(' ')[0] for s in allSpecs]))

spec_taxids = ncbi.get_name_translator(allSpecs)

gen_taxids = ncbi.get_name_translator(genera)

# for s in allSpecs:
#     if s not in spec_taxids.keys():
#         print(s)


# for g in genera:
#     if g not in gen_taxids.keys():
#         print(g)
#         for s in allSpecs:
#             if s.split(' ')[0] == g:
#                 print('\t' + s)
                
# Check non-viral genera for possible one-off spelling errors
name_changes = {} # maps from corrected species name(most specific level provided) to original name

# Kinetoplastids -> Kinetoplastida
# ncbi.get_name_translator(['Kinetoplastida'])
name_changes['Kinetoplastida'] = 'Kinetoplastids'

# Actinogyra muehlenbergii -> Umbilicaria muehlenbergii
## Actinogyra is deprecated
# ncbi.get_name_translator(['Umbilicaria muehlenbergii'])
# ncbi.get_name_translator(['Umbilicaria'])
name_changes['Umbilicaria muehlenbergii'] = 'Actinogyra muehlenbergii'

# Arecastrum romanzoffianum -> Syagrus romanzoffiana
## Queen palm
# ncbi.get_name_translator(['Syagrus'])
# ncbi.get_name_translator(['Syagrus romanzoffiana'])
name_changes['Syagrus romanzoffiana'] = 'Arecastrum romanzoffianum'

# Columbia livia -> Columba livia
## Rock dove (pigeon)
# ncbi.get_name_translator(['Columba'])
# ncbi.get_name_translator(['Columba livia'])
name_changes['Columba livia'] = 'Columbia livia'


genera.extend([s.split(' ')[0] for s in name_changes.keys() if s.split(' ')[0] not in genera])

gen_taxids = ncbi.get_name_translator(genera)

for k,v in name_changes.items(): # Duplicate entries for missnamed species to still map to correct taxid
    gen_taxids[v.split(' ')[0]] = gen_taxids[k.split(' ')[0]]
    
# gen_taxids['Columba']
# gen_taxids['Columbia']

# lineage = ncbi.get_rank(ncbi.get_lineage(gen_taxids['Columbia'][0]))

# lin_names = ncbi.get_taxid_translator(list(lineage.keys()))

# for k,v in lineage.items():
#     if v == 'superkingdom':
#         print(lin_names[k])

# for k,v in lineage.items():
#     if v == 'kingdom':
#         print(lin_names[k])
        
# for k,v in lineage.items():
#     if v == 'phylum':
#         print(lin_names[k])


# Get taxonomic classifiers for each taxid
superkingdoms = collections.defaultdict(list)
kingdoms = collections.defaultdict(list)
phyla = collections.defaultdict(list)
classes = collections.defaultdict(list)
orders = collections.defaultdict(list)
families = collections.defaultdict(list)


for genus in gen_taxids.keys():
    lineage = ncbi.get_rank(ncbi.get_lineage(gen_taxids[genus][0])) # Get lineage for current genus and translate to ranks
    lin_names = ncbi.get_taxid_translator(list(lineage.keys())) # Get names for each taxid in lineage
    
    for k,v in lineage.items():
        if v == 'superkingdom':
            superkingdoms[genus].append(lin_names[k])
        if v == 'kingdom':
            kingdoms[genus].append(lin_names[k])
        if v == 'phylum':
            phyla[genus].append(lin_names[k])
        if v == 'class':
            classes[genus].append(lin_names[k])
        if v == 'order':
            orders[genus].append(lin_names[k])
        if v == 'family':
            families[genus].append(lin_names[k])
        
    
uni_SKs = list(superkingdoms.values())
uni_Ks = list(kingdoms.values())
uni_Ps = list(phyla.values())
uni_Cs = list(classes.values())
uni_Os = list(orders.values())
uni_Fs = list(families.values())
 
####

sum([1 for sk in uni_SKs if len(sk) == 1]) == len(uni_SKs)
len(uni_SKs)
sum([1 for k in uni_Ks if len(k) == 1]) == len(uni_Ks)
len(uni_Ks)
sum([1 for p in uni_Ps if len(p) == 1]) == len(uni_Ps)
len(uni_Ps)
sum([1 for t in uni_Cs if len(t) == 1]) == len(uni_Cs)
len(uni_Cs)
sum([1 for t in uni_Os if len(t) == 1]) == len(uni_Os)
len(uni_Os)
sum([1 for t in uni_Fs if len(t) == 1]) == len(uni_Fs)
len(uni_Fs)

## No genera taxids map to multiple taxonomic classifiers of the same rank
## All 826 genera have a superkingdom
# 233 are missing kingdoms
# 4 are missing phyla
# 3 are missing classes
# 9 are missing orders
# 18 are missing families

# Only eukaryotes are assigned to kingdoms, add Archaea, Bateria, & Virus superkingdoms as separate kingdoms 
## (virus labels added to superkingdoms (& now kingdoms) in the encoding script)
for k,v in superkingdoms.items():
    if k not in kingdoms and v[0] in ['Bacteria', 'Archaea']:
        kingdoms[k].extend(v)

uni_Ks = list(kingdoms.values())
sum([1 for k in uni_Ks if len(k) == 1]) == len(uni_Ks)
len(uni_Ks)
# Now only 14 missing kingdoms
####

uni_SKs = [sk[0] for sk in uni_SKs]
uni_Ks = [k[0] for k in uni_Ks]
uni_Ps = [p[0] for p in uni_Ps]
uni_Cs = [t[0] for t in uni_Cs]
uni_Os = [t[0] for t in uni_Os]
uni_Fs = [t[0] for t in uni_Fs]

collections.Counter(uni_SKs) # 3 superkingdoms(domains)
collections.Counter(uni_Ks) # 5 kingdoms
len(set(uni_Ps)) # 28 phyla
len(set(uni_Cs)) # 84 classes
len(set(uni_Os)) # 207 orders
len(set(uni_Fs)) # 408 families


# Save genus to taxanomic rank dictionaries

taxa_relations = [superkingdoms, kingdoms, phyla, classes, orders, families]

# Since mappings are 1:1, remove list format in values
for dic in taxa_relations:
    for g,t in dic.items():
        dic[g] = t[0]

with open('/Users/dmattox/cbk/sugarTrees/pickles/genera2ranks.pickle', 'wb') as outFH:
    pickle.dump(taxa_relations, outFH, protocol=pickle.HIGHEST_PROTOCOL)



