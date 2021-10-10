#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:02:43 2021

@author: dmattox
"""

import re, collections, pickle
import numpy as np

import cbk.sugarTrees.bin.sugarTrees as SugarTrees
import cbk.sugarTrees.bin.getStats as SugarStats

homeDir = '/Users/dmattox/cbk/sugarTrees/'

cfgFile = homeDir+'data/Clean_Lectin_intensities_67uni.txt'

cfgGlys = []
cfgGly2cfgID = collections.defaultdict(list)

with open(cfgFile, 'r') as inFH: # Read in sugar base csv file
    for line in inFH:
        line = line.strip()
        line = line.split('\t')
        if line[0] != 'glyID':
            cfgGlys.append(line[1].strip())
            cfgGly2cfgID[line[1].strip()].append(int(line[0]))


# monosaccGlys = [g.split('-')[0] for g in cfgGlys if len(g.split('-')) == 2] # monosaccharides as individual glycans
# list(set([m[:-1] for m in monosaccGlys]))



print(len(cfgGlys)) # 610 glycans

# Check out linker variety and frequency
allLinkers = [re.findall('-Sp[0-9]*', gly, flags = re.IGNORECASE) for gly in cfgGlys]
for i,l in enumerate(allLinkers):
    if len(l) != 1:
        print(i, l, cfgGlys[i])

allLinkers = [l for lst in allLinkers for l in lst]
collections.Counter(allLinkers)

noAno_linkers = [re.search('-Sp[0-9]*', gly, flags = re.IGNORECASE)[0] for gly in cfgGlys if re.sub(r'-Sp[0-9]*','', gly, flags = re.IGNORECASE)[-1] not in ['a', 'b']]
collections.Counter(noAno_linkers) # some linkers have defined anomeric conformations, but none of these linkers where anomeric conf. is not specified
## https://doi.org/10.1093/glycob/cwt083


for gly in cfgGlys:
    linkless = re.sub(r'-Sp[0-9]*','', gly, flags = re.IGNORECASE)
    cfgGly2cfgID[linkless].extend(cfgGly2cfgID[gly])
    
cfgGlys = list(set([re.sub('-Sp[0-9]*','', gly, flags = re.IGNORECASE) for gly in cfgGlys])) # trim the linker info
print(len(cfgGlys)) # 550 unique glycans


# # Leave out glycans composed of single monosaccharides?
# cfgGlys = [gly for gly in cfgGlys if len(gly.split('-')) >= 2]
# print(len(cfgGlys)) # 529 with anomeric info

for k in list(cfgGly2cfgID.keys()):
    if k not in cfgGlys:
        del cfgGly2cfgID[k]
        
len(cfgGly2cfgID)

# Double check that links were removed and glycans have been process correctly
collections.Counter([g[-1] for g in cfgGlys])
[g for g in cfgGlys if g[-1] == 'n']
[g for g in cfgGlys if g[-1] == 'l']
[g for g in cfgGlys if g[-1] == '1'] # two glycans have a carbon linkage 

# Check out frequency of defined anomeric state of final sugar
sum([1 for gly in cfgGlys if gly[-1] in ['a', 'b', '1']])

no_ano = [gly[:-1] if gly[-1] in ['a', 'b'] else gly[:-2] if gly[-1] == '1' else gly for gly in cfgGlys]
len(no_ano)

ano_cnts = []
for gly in cfgGlys:
    if gly[-1] not in ['a', 'b']:
        inds = [i for i,g in enumerate(no_ano) if g == gly]
        if len(inds) > 1:
            ano_cnts.append(len(inds))
print(len(ano_cnts))

# Reformat anomeric info of root saccharide into empty links
for gly in cfgGlys:
    newGly = re.sub(r'([a|b]$)', r'\g<0>0-0', gly)
    newGly = re.sub(r'([a|b]1$)', r'\g<0>-0', newGly) # edge case where the carbon connected to the linker on the root sugar is specified
    cfgGly2cfgID[newGly].extend(cfgGly2cfgID[gly])
    
len(cfgGly2cfgID)
    
cfgGlys = [re.sub(r'([a|b]$)', r'\g<0>0-0', gly) for gly in cfgGlys]
cfgGlys = [re.sub(r'([a|b]1$)', r'\g<0>-0', gly) for gly in cfgGlys]

for k in list(cfgGly2cfgID.keys()):
    if k not in cfgGlys:
        del cfgGly2cfgID[k]

len(cfgGly2cfgID)


# get monosaccharide vocabulary
cfg_monosaccs = []
for gly in cfgGlys:
    sacc_list = re.split('[a|b][0-9]-[0-9]', gly) # split monosaccharides by linkage (removes linkage)
    
    # drop empty string leftover from pseudolink of root sacc
    sacc_list = [sacc for sacc in sacc_list if sacc != '']
    
    sacc_list = [sacc[1:] if sacc[0] == ')' else sacc for sacc in sacc_list] # Trim closing parenthesis leftover from branches if at beginning of sugar name
    sacc_list = [sacc[1:] if ((sacc[0] == '(') and (len(re.findall('\(', sacc)) != len(re.findall('\)', sacc)))) else sacc for sacc in sacc_list] # Trim opening parenthesis leftover from branches if at beginning of sugar name and isn't matched by a closing parenthesis
    
    sacc_list = [sacc.strip() for sacc in sacc_list] # strip whitespace
    sacc_list = [sacc[:-1] if sacc[-1] == '-' else sacc for sacc in sacc_list] # Strip extraneuos linkage info (should only be one instance)
    
    if 'GlcNAcb1' in sacc_list:
        print(gly)
        print(cfgGly2cfgID[gly])

    cfg_monosaccs.extend(sacc_list)

cfg_mcnts = collections.Counter(cfg_monosaccs)
cfg_mcnts.most_common()

cfg_monosaccs = list(set(cfg_monosaccs))
print(len(cfg_monosaccs)) # 33 unique sugars in these glycans

newSugars = [gly for gly in cfg_monosaccs if gly not in SugarStats.commonMs]
print(len(newSugars)) # 24 of these sugars are not seen in the more common sugarbase monosaccharides

len([gly for gly in cfg_monosaccs if gly not in SugarStats.monosaccLst]) # Also not found in any of the monosaccharides from SugarBase, likely a formatting issue

# Hierarchy of modifications (https://doi.org/10.1101/2020.04.08.031948)
##  NAc > OAc > NGc > OGc > NS > OS > NP > OP > NAm > OAm > NBut > OBut > NProp > OProp > NMe > OMe > CMe > NFo > OFo > OPPEtn > OPEtn > OEtn > A > N 
##  > SH > OPCho > OPyr > OVac > OPam > OEtg > OFer > OSin > OAep > OCoum > ODco > OLau > OSte > OOle > OBz > OCin > OAch > OMal > OMar > OOrn > rest.

# reconcile misformatted vocabulary between sugarbase and CFG microarrays
## NOTE: Original reconciliation performed with monosaccharides from monosaccharide glycans as well
cfg2SB_monosaccs = dict.fromkeys(newSugars)

SugarStats.commonMs

##
[k for k,v in cfg2SB_monosaccs.items() if v is None] # Run to see list of monosacchrides left to process
##

# KDN --> Kdn
'Kdn' in SugarStats.commonMs
cfg2SB_monosaccs['KDN'] = 'Kdn'


# Neu5Gc --> NeuNGc
[g for g in SugarStats.commonMs if g[:3] == 'Neu']
cfg2SB_monosaccs['Neu5Gc'] = 'NeuNGc'

# Neu5Ac --> NeuNAc
cfg2SB_monosaccs['Neu5Ac'] = 'NeuNAc'

# Neu5,9Ac2 --> NeuNAcOAc (5-N-acetyl-9-O-acetyl neuraminic acid)
sorted([e for e in SugarStats.commonMs if e[:3] == 'Neu'])
'NeuNAcOAc' in SugarStats.commonMs
cfg2SB_monosaccs['Neu5,9Ac2'] = 'NeuNAcOAc'


# (6S)Glc --> GlcOS
'GlcOS' in SugarStats.commonMs
sorted([e for e in SugarStats.commonMs if e[:3] == 'Glc'])
cfg2SB_monosaccs['(6S)Glc'] = 'GlcOS'

# (6P)Glc --> GlcOP
'GlcOP' in SugarStats.commonMs
cfg2SB_monosaccs['(6P)Glc'] = 'GlcOP'

# (3S)GlcA --> GlcOSA
'GlcOSA' in SugarStats.commonMs
cfg2SB_monosaccs['(3S)GlcA'] = 'GlcOSA'


# (6S)GlcNAc --> GlcNAcOS
'GlcNAcOS' in SugarStats.commonMs
cfg2SB_monosaccs['(6S)GlcNAc'] = 'GlcNAcOS'

# GlcNA --> GlcAN
'GlcAN' in SugarStats.commonMs
cfg2SB_monosaccs['GlcNA'] = 'GlcAN'

# GlcN(Gc) --> GlcNGc
'GlcNGc' in SugarStats.commonMs
'GlcNGc' in SugarStats.monosaccLst # GlcNGc falls into less common glycan category
########
del cfg2SB_monosaccs['GlcN(Gc)']
[gly for gly in cfgGlys if gly.find('GlcN(Gc)') != -1]
cfgGlys.remove('GlcN(Gc)b0-0')
del cfgGly2cfgID['GlcN(Gc)b0-0']
print(len(cfgGly2cfgID), len(cfgGlys))
########

# (3S)GlcNAc --> GlcNAcOS
'GlcNAcOS' in SugarStats.commonMs
cfg2SB_monosaccs['(3S)GlcNAc'] = 'GlcNAcOS'

# (6P)GlcNAc --> GlcNAcOP
'GlcNAcOP' in SugarStats.commonMs
cfg2SB_monosaccs['(6P)GlcNAc'] = 'GlcNAcOS'


# (6P)Man --> ManOP
sorted([e for e in SugarStats.commonMs if e[:3] == 'Man'])
'ManOP' in SugarStats.commonMs
cfg2SB_monosaccs['(6P)Man'] = 'ManOP'


# (6P)Gal  --> GalOP
'GalOP' in SugarStats.commonMs
cfg2SB_monosaccs['(6P)Gal'] = 'GalOP'


# (3/4/6S)GalNAc --> GalNAcOS
sorted([e for e in SugarStats.commonMs if e[:3] == 'Gal']) #sugarbase has no distinctions for positions or degree of sulfation, collapse all to "OS"
sorted([e for e in SugarStats.monosaccLst if e[:3] == 'Gal'])
'GalNAcOS' in SugarStats.commonMs
cfg2SB_monosaccs['(6S)GalNAc'] = 'GalNAcOS'
cfg2SB_monosaccs['(4S)GalNAc'] = 'GalNAcOS'
cfg2SB_monosaccs['(3S)GalNAc'] = 'GalNAcOS'
cfg2SB_monosaccs['(6S)(4S)GalNAc'] = 'GalNAcOS'

# all remaining sulfated and disulfated galctoses --> GalOS
'GalOS' in SugarStats.commonMs
cfg2SB_monosaccs['(3S)Gal'] = 'GalOS'
cfg2SB_monosaccs['(4S)Gal'] = 'GalOS'
cfg2SB_monosaccs['(6S)Gal'] = 'GalOS'

cfg2SB_monosaccs['4S(3S)Gal'] = 'GalOS'
cfg2SB_monosaccs['6S(3S)Gal'] = 'GalOS'
cfg2SB_monosaccs['(6S)(4S)Gal'] = 'GalOS'



cfg2SB_monosaccs


# Reformat branches
len(cfgGlys)
clean_cfgGlys = [] # hold the final reformatting of CFG glycans

for gly in cfgGlys: # Making sure there aren't parenthesis outside of the misformatted monosaccharides and the branches
    sacc_list = re.split('[a|b][0-9]-[0-9]', gly)
    sacc_list = [sacc for sacc in sacc_list if sacc != ''] # drop empty string leftover from pseudolink of root sacc
    sacc_list = [sacc[1:] if sacc[0] == ')' else sacc for sacc in sacc_list] # Trim closing parenthesis leftover from branches if at beginning of sugar name
    sacc_list = [sacc[1:] if ((sacc[0] == '(') and (len(re.findall('\(', sacc)) != len(re.findall('\)', sacc)))) else sacc for sacc in sacc_list] # Trim opening parenthesis leftover from branches if at beginning of sugar name and isn't matched by a closing parenthesis

    for s in sacc_list:
        if s.find('(') != -1 and s not in cfg2SB_monosaccs.keys():
            print(gly)
            print(s)
            print()
            
# Replaces parentheses around branches with square brackets    
clean_cfgGlys = [re.sub('\(', '[', gly) for gly in cfgGlys]
clean_cfgGlys = [re.sub('\)', ']', gly) for gly in clean_cfgGlys]

# Duplicated CFG saccharides with parenthesis in dictionary to represent brackets
newKs = [(k,v) for k,v in cfg2SB_monosaccs.items() if k.find('(') != -1]
for k,v in newKs:
    k = re.sub('\(', '[', k)
    k = re.sub('\)', ']', k)
    cfg2SB_monosaccs[k] = v
    
# Place parentheses around links
clean_cfgGlys = [re.sub(r'([a|b][0-9]-[0-9])', r'(\1)', gly) for gly in clean_cfgGlys]


# Clean extranous hyphens, whitespaces in sugar names, and swap out CFG monosaccharides for SB ones
for i,gly in enumerate(clean_cfgGlys):
    links = re.findall(r'\([a|b][0-9]-[0-9]\)', gly) # get all linkages in original order
    sacc_list = re.split(r'\([a|b][0-9]-[0-9]\)', gly) # Get all sugars and branch brackets
    sacc_list = [sacc for sacc in sacc_list if sacc != ''] # drop empty string leftover from pseudolink of root sacc
    
    # Pull out closing brackets leftover from branches
    sacc_list = [sacc[0] + '!' + sacc[1:] if (sacc[0] == ']') else sacc for sacc in sacc_list] # Add exclaimation marks after unmatched closing brackets in a monosaccharide
    sacc_list = '!'.join(sacc_list).split('!')
    # Pull out opening brackets leftover from branches
    sacc_list = [sacc[0] + '!' + sacc[1:] if ((sacc[0] == '[') and (len(re.findall('\[', sacc)) != len(re.findall('\]', sacc)))) else sacc for sacc in sacc_list] # Add exclaimation marks after unmatched opening brackets in a monosaccharide
    sacc_list = '!'.join(sacc_list).split('!')
    
    sacc_list = [sacc.strip() for sacc in sacc_list] # strip whitespace
    sacc_list = [sacc[:-1] if sacc[-1] == '-' else sacc for sacc in sacc_list] # Strip extraneuos linkage info (should only be one instance, see next line)
    ## Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man-(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc
    
    # confirm all sugars are correctly formatted or accounted for
    if any([True if (s not in SugarStats.commonMs and s not in cfg2SB_monosaccs.keys() and s not in ['[', ']']) else False for s in sacc_list ]):
        print(sacc_list)
    
    # swap CFG formatted sugars for SB formatting
    sacc_list = [cfg2SB_monosaccs[s] if s in cfg2SB_monosaccs.keys() else s for s in sacc_list]

    # Reconstruct glycan
    newGly = []
    
    while links or sacc_list:
        if not newGly:
            newGly.append(sacc_list.pop(0))
            
        if newGly[-1][-1] != ')' and newGly[-1] not in ['[',']']:
            newGly.append(links.pop(0))
        else:
            newGly.append(sacc_list.pop(0))
        
    clean_cfgGlys[i] = ''.join(newGly)

# Check reduced number of unique glycans after cleaning
[g for i,g in enumerate(clean_cfgGlys) if g in cfgGly2cfgID.keys()]
len([g for i,g in enumerate(cfgGlys) if g in cfgGly2cfgID.keys()])
[g for i,g in enumerate(cfgGlys) if g not in cfgGly2cfgID.keys()]
    
len(clean_cfgGlys)
len(set(clean_cfgGlys)) # Lack of detail in OS 
len(cfgGlys)

dupGlys = collections.defaultdict(list)
for i,g in enumerate(clean_cfgGlys):
    if clean_cfgGlys.count(g) > 1:
        dupGlys[g].append(cfgGlys[i])
len(dupGlys.keys())

cnt = 0
for k,v in dupGlys.items():
    print(k)
    for g in v:
        cnt += 1
        print('\t' + g)
    print()
print(str(cnt)+' CFG glycans collapse to '+str(len(dupGlys.keys()))+' cleaned SB-formatted glycans' )

# Update dictionary
for i,g in enumerate(cfgGlys):    
    cfgGly2cfgID[clean_cfgGlys[i]].extend(cfgGly2cfgID[g])
    del cfgGly2cfgID[g]

len(clean_cfgGlys)
len(cfgGly2cfgID)

clean_cfgGlys = list(set(clean_cfgGlys)) # 
print(len(clean_cfgGlys))

for k,v in cfgGly2cfgID.items():
    cfgGly2cfgID[k] = sorted(list(set(v)))


cfg_monosaccs = []
for gly in clean_cfgGlys:
    sacc_list = re.split('[a|b][0-9]-[0-9]', gly) # split monosaccharides by linkage (removes linkage)
    
    # trim link parentheses and branch brackets
    sacc_list = [re.sub(r'[(|)|\[|\]]', '', sacc) for sacc in sacc_list]
    
    # drop empty strings leftover from branches and links
    sacc_list = [sacc for sacc in sacc_list if sacc != '']

    cfg_monosaccs.extend(sacc_list)

cfg_mcnts = collections.Counter(cfg_monosaccs)
cfg_mcnts.most_common()

len(set(cfg_monosaccs))


len([gly for gly in clean_cfgGlys if gly in list(cfgGly2cfgID.keys())])


cfgGlys2sbIDs = dict.fromkeys(clean_cfgGlys)

for gly in clean_cfgGlys:
    for sbEntry in SugarStats.okayEntries:
        if gly == sbEntry.iupac:
            cfgGlys2sbIDs[gly] = sbEntry.sbID

len([k for k,v in cfgGlys2sbIDs.items() if v is not None]) # 39 of 408 CFG glycans accounted for

# Some might be missing due to shuffled branches, others missing becasue of added anomeric info
rootLinks = {}
for gly in clean_cfgGlys:
    if  re.search(r'\([a|b][0-9]-[0-9]\)$', gly):
        rootLinks['_'.join([str(id) for id in cfgGly2cfgID[gly]])] = re.search(r'\([a|b][0-9]-[0-9]\)$', gly)[0]
    

cfgGlyEncodings = {}
for gly in clean_cfgGlys:
    cfgGlyEncodings['_'.join([str(id) for id in cfgGly2cfgID[gly]])] = SugarTrees.SugarBase(*['_'.join([str(id) for id in cfgGly2cfgID[gly]]), re.sub(r'\([a|b][0-9]-[0-9]\)$', '', gly), 'None', '', 'Unknown'])



for k,gly in cfgGlyEncodings.items():
    try:
        gly.buildTree()
        gly.drawSbwyMap()
        if k in rootLinks:
            link = SugarTrees.getLink(rootLinks[k])[1]
            gly.tree[(1,0)].anomeric = link[0] if link[0] in ['a', 'b'] else 'u' # update anomeric conformation of root glycan (reducing end) if stored in pseudolink
            link = link[1:] if link[0] in ['a', 'b', 'u'] else link
            if link != '0-0': # '0-0' is the default link used to build the root node, if it's the same there is nothing else to update
                carbonLinks = link.split('-')
                if carbonLinks[0] != '0':
                    del gly.tree[(1,0)].Clinks['0'] # drop carbon linkage info from default rootnode info
                    gly.tree[(1,0)].Clinks[carbonLinks[0]] = 2 # add carbon linkage info from provided pseudolink
    except:
        print('Warning\n\t CFG glycan ' + k, gly.iupac, 'failed to process')
        gly.tree = {}

del cfgGlyEncodings['509']

print(len(cfgGlyEncodings))

#############################
# Get encodings for Bowen
#############################

# Make sure everything is within ranges observed for set of SugarBase glycans

# All monosacchardies are shared and accounted for
monosaccLst = set()
for sugar in cfgGlyEncodings.values():
    monosaccLst.update([s.base for s in sugar.tree.values()])
monosaccLst = list(monosaccLst)
len(monosaccLst)
[m for m in monosaccLst if m not in SugarStats.monoDict.keys()]


# Max branches in CFG glycans is 6 (11 from SB glycans)
maxB = [max(sugar.tree[(0,0)].subway) for sugar in cfgGlyEncodings.values()]

cntB = collections.Counter(maxB)
cntB.most_common()

maxB = max(maxB)
if maxB > SugarStats.maxB:
    raise Warning('Highest observed degree of branching is greater than observed in all glycans used for pre-processing, encodings will not be entirely consistent')
else:
    maxB = SugarStats.maxB


# Max carbon linkage in CFG glycans is 8 (9 from SB glycans)
maxC = []
for sugar in cfgGlyEncodings.values():
    try:
        maxC.append(int(max([max(s.Clinks.keys()) for s in sugar.tree.values()])))
        if int(max([max(s.Clinks.keys()) for s in sugar.tree.values()])) == 0:
            print(sugar.iupac) # Should only be monosaccharides
    except:
        print(sugar)
        
cntC = collections.Counter(maxC)
cntC.most_common()
maxC = max(maxC)

if maxC > SugarStats.maxC:
    raise Warning('Highest observed carbon number involved in a glycosidic bond is greater than observed in all glycans used for pre-processing, encodings will not be entirely consistent')
else:
    maxC = SugarStats.maxC




cfgGlycanEncodings_535_apr12 = {}

# colNames = ['sugar'] + ['anomeric'] + ['C_links'] + ['B_lines'] + ['parLink', 'sDepth', 'sInd']

topColVals = [0]*7

for k,sugar in cfgGlyEncodings.items():
    
    sugar.buildEncoding(maxB = SugarStats.maxB, maxC = SugarStats.maxC, monoDict = SugarStats.monoDict)
    
    for i in range(len(topColVals)):
        if max(sugar.encoding[:,i]) > topColVals[i]:
            topColVals[i] = max(sugar.encoding[:,i])
        
    li = SugarStats.linkDict[sugar.link] # Hold link code
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
    im = SugarStats.immDict[sugar.immunogenic]
    
    cfgGlycanEncodings_535_apr12[k] = [sugar.encoding, li, sp, im]


print(len(cfgGlycanEncodings_535_apr12))

print(topColVals)

# sugar
# sugar.print()
# sugar.buildEncoding(maxB = maxB, maxC = maxC, monoDict = SugarStats.monoDict)
# sugar.encoding

# np.array_equal(sugar.encoding, cfgGlycanEncodings_535_apr12['102_103'][0])

# Serialize encodings for Bowen
encodeFile = homeDir+'pickles/CFG_glycans_535_apr12.pickle'

with open(encodeFile, 'wb') as outFH:
    pickle.dump(cfgGlycanEncodings_535_apr12, outFH, protocol=pickle.HIGHEST_PROTOCOL)




# Save out processed SugarTree objects too
encodeFile = homeDir+'pickles/CFG_glycans_535_sugarTrees.pickle'

with open(encodeFile, 'wb') as outFH:
    pickle.dump(cfgGlyEncodings, outFH, protocol=pickle.HIGHEST_PROTOCOL)

##############################
# Read in and process embeddings into .tsv file
##############################

embFile = homeDir + 'pickles/CFG_glycans_535_feb9_emb.pickle'

with open(embFile, 'rb') as inFH:
    cfgGlyEmbeds = pickle.load(inFH)

# len(cfgGlyEmbeds)
# cfgGlyEmbeds['190'].shape

# cfgGlyEncodings['190'].iupac
# cfgGlyEncodings['190'].tree

# cfgGlyEmbeds['190'][0,:,0]

colnames = ['glyID'] + [i + j for i, j in zip(['d']*384, [str(i) for i in range(384)])]

embCSV = homeDir + 'CFG_glycans_535_feb9_emb.csv'

with open(embCSV, 'w') as outFH:
    outFH.write(','.join(colnames) + '\n')
    for k,gly in cfgGlyEmbeds.items():
        emb = gly[0,:,0]
        outFH.write(','.join([k] + [str(d) for d in emb]) + '\n')



































