#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 12:26:18 2020

@author: dmattox
"""

# import numpy as np
import collections
import numpy as np

from anytree import NodeMixin, RenderTree


def binary2decimal(binary):
    # Accepts binary number as a string, returns decimal number as integer
    out = 0
    binary = binary[::-1] # Reverse direction of binary string
    
    for i, b in enumerate(binary):
        out += int(b) * (2 ** i)
        
    return out
        
def decimal2binary(decimal):
    # Accepts decimal number as integer, returns binary number as a string
    out = ''
    
    while decimal > 0:
        out += str(decimal % 2) # bit equal to remainder of dividing by 2
        decimal //= 2 # update value to iteger quotient
        
    return out[::-1] # Reverse order of binary string and return it

def validParenth(gly):
    # Returns true if all parentheses and brackets are paired and closed
    stack = []
    comp = {')': '(',
           ']': '['}
    for p in gly:
        if p in ['(',')','[',']','{','}']:
            if p in comp.values():
                stack.insert(0,p)
            else:
                if not stack:
                    return False
                elif stack.pop(0) != comp[p]:
                    return False
    if stack:
        return False
    else:
        return True

def getSug(iupacGly):
    mono = ''
    for l in iupacGly[::-1]: # loop through string backwards and pop off first complete monosaccharide
        if l not in ['[', ']', '(', ')']:
            mono += l
        else:
            break
    if mono == '':
        raise SyntaxError('ERROR: Linkage punctuation found before any sugars, returning unchanged iupac glycan string')
    mono = mono[::-1] # correct order for sugar iupac code
    return(iupacGly[:(-1*len(mono))], mono)

def getLink(iupacGly):
    link = ''
    for l in iupacGly[::-1]: # loop through string backwards and pop off the link defined within the last closed parenthesises
        if l != '(':
            link += l
        else:
            link += l
            break
    link = link[::-1] # correct order for sugar iupac code
    return(iupacGly[:(-1*len(link))], link[1:-1]) # Return iupac string without link, as well as the link without the parentheses

def getBranch(iupacGly):
    # Pop an entire  bracketed branch from the iupac string
    branch = ''
    nested = []
    for l in iupacGly[::-1]:
        if l == '[':
            if nested == []:
                branch += l
                break
            else:
                branch += l
                nested.pop(0)
        elif l == ']' and branch != '': # If it hits a nested branch (BESIDES THE FIRST BRACKET)
            branch += l
            nested.append(l)
        else:
            branch += l
    branch = branch[::-1] # Reverse back to correct order
    return(iupacGly[:(-1*len(branch))], branch[1:-1]) # Return iupac string without branch, as well as the branch without the parentheses

# def decode(self, maxB, maxC, monoDict, anoDict = {'u': 0, 'a': 1, 'b': 2}):
#     """
    

#     Parameters (parameters used to generate encoding)
#     ----------
#     maxB : int
#         Max value contained in a subway line in the root nodes of all glycans being considered (maximum number of observed branches).
#     maxC : int
#         Highest carbon number observed to participate in a glycosidic bond (from all glycans being considered).
#     monoDict : dict
#         Dictionary linking all observed monosaccharides to a corresponding integer (integers based on monosaccharide frequency rank).
#     anoDict : dict, optional
#         Dictionary encoding anomeric conformation information as integers.  The default is {'u': 0, 'a': 1, 'b': 2}.

#     Returns
#     -------
#     IUPAC string from encoding

#     """
#     pass


class GlyNode(NodeMixin):
    # Class to hold nodes (rows) of tree representations
    def __init__(self, base, ind, link = 'u0-0', parent = None):
        super(GlyNode, self).__init__()
        self.base = base
        # self.depth = depth # Depth is an existing property of parent class
        self.ind = ind
        
        link = link.replace('(','') # clean extra parentheses
        link = link.replace(')','')
        
        self.anomeric = link[0] if link[0] in ['u','a','b'] else 'u' # state of anomeric carbon if indicated, otherwise assume unknown
        self.parent = parent # parent GlyNode

        self.name = '%s:%d*%d' % (self.base, self.depth, self.ind)
        
        link = link[1:] if link[0] in ['u','a','b'] else link[:] # drop anomeric charcter if present
        link = link.split('-') # Split linkage into child Carbon connection [0] and parent carbon connection [1]
        
        self.Clinks = collections.Counter() # Holds the attachemnts at different carbon positions (based on formal carbon numbering)
        
        if self.parent is None:
            self.parentLink = 0
        else:
            self.Clinks[link[0]] = 2
            self.parent.Clinks[link[1]] = 1
            self.parentLink = link[1]
            # self.parent.children.append(self) # Children handled by anytree parent class
 
        
        # self.children = [] # handled by anytree parent class
        self.subway = [] # Holds the indices of the subway lines that pass through the node
        
    def __repr__(self):
        return self.name
    def __str__(self):
        return self.name
    
    def treePrint(self):
        for pre, fill, node in RenderTree(self):
            print("%s%s" % (pre, node.name))
    
    def drawSubway(self, sbwyStp):
        # Called by SugarBase.drawSbwyMap
        
        self.subway.append(sbwyStp) # Add subway line for the terminal node
        if self.parent is not None:
            self.parent.drawSubway(sbwyStp) # Recursively pass the subway line up the tree towards the root node, adding the information to each node along the way
    
class SugarBase:
    # Class to hold entires in SugarBase v2 database
    def __init__(self, sbID, iupac, link, species, immunogenic):
        self.sbID = sbID
        self.iupac = iupac
        self.link = link
        self.species = [s.strip() for s in species.split(',')]
        self.taxonomy = []
        self.immunogenic = immunogenic
        
        self.id = int(sbID.replace('SBID', ''))
        
        self.tree = {}
        # self.buildTree()
        
        self.encoding = None
    
    def  __repr__(self):
        return self.sbID
    
    def __str__(self):
        return self.sbID
    
    def __int__(self):
        return(self.id)
    
    def print(self):
        print('SugarBase ID:\n\t', self.id)
        print('N/O Linked:\n\t', self.link)
        if self.species[0] == '':
            print('Species of origin:\n\tUnknown')
        elif len(self.species) <= 2:
            print('Species of origin: [', len(self.species),']\n\t', self.species)
        else:
            print('Species of origin: [', len(self.species),']\n\t', self.species[:3], '... (first 3, try print(SugarBase.species) to see the rest )')
        print('Immunogenicity:\n\t', self.immunogenic)
        print('IUPAC glycan:\n\t', self.iupac)
    
    def treePrint(self):
        if self.tree == {}:
            print('ERROR: sugar tree not yet constructed')
        else:
            self.tree[(0,0)].treePrint()
    
    def treeDepthCnt(self, depth):
        # Returns the number of nodes in the tree at the given depth
        cnt = 0
        for k in self.tree.keys():
            if k[0] == depth:
                cnt += 1
        return cnt
    
    def buildTree(self):
        if self.tree != {}:
            print('WARNING: Tree already constructed, not rebuilding it')
        elif (validParenth(self.iupac) == False):
            raise SyntaxError('Unmatched parentheses or brackets detected in supplied IUPAC glyan string')
        else:
            gly = self.iupac

            # add start token
            par = (0,0)
            self.tree[par] = GlyNode(base = 'START', ind = par[1])
            
            
            # Process the root node
            gly, base = getSug(gly)
            chi = (1,0)
            self.tree[chi] = GlyNode(base = base, ind = par[1], parent = self.tree[par])
            par = chi
            
            if gly: # if glycan is a monosaccharide, sets the queue to empty to avoid the while loop
                branchQueue = [[gly,par]]
            else:
                branchQueue = []
            
            while branchQueue:
                
                if branchQueue[0][0][-1] != ')' and branchQueue[0][0][-1] != ']':
                    print('ERROR: no linkage or branch found for glycan ', self.sbID)
                    break
            
                if branchQueue[0][0][-1] == ']':
                    par = branchQueue[0][1]
                    childLst = [] # Branching, at least 2 children from current parent node
                    while branchQueue[0][0][-1] == ']':
                        branchQueue[0][0], branch = getBranch(branchQueue[0][0])
                        childLst.append(branch)
                        
                    childLst.append(branchQueue[0][0]) # add "main" branch to the list of branches as well 
                    branchQueue.pop(0) # and remove it from the original queue
                    
                    childLst.sort(key = lambda x: int(x[-2]), reverse = True) # sort all of the branches from the parent node by descending parentlink carbon number
                    
                    for branch in childLst:
                        branchQueue.insert(0,[branch, par]) # Add braches to branch queue such that the lower numbered branches are on the top of the queue
                    chi = par # Since no monosacchairdes are removed, set chi to par to preserve true parent
            
                if branchQueue[0][0][-1] == ')':
                    
                    par = branchQueue[0][1]
                    chi = (par[0]+1, self.treeDepthCnt(par[0]+1)) # depth & index of child
                    
                    branchQueue[0][0], link = getLink(branchQueue[0][0])
                    branchQueue[0][0], base = getSug(branchQueue[0][0])
                
                    self.tree[chi] = GlyNode(base, ind=chi[1], link=link, parent = self.tree[par])
            
                
                
                if branchQueue[0][0] == '':
                    branchQueue.pop(0) # If a branch has been fully processed, remove it from the queue
                else:
                    branchQueue[0][1] = chi # otherwise, update the parent for the remainder of the branch
            # Add stop tokens to terminal monosaccharides
            termNodes = []
            for k,v in self.tree.items():
                if v.children == ():
                    termNodes.append(k)
            termNodes.sort(key= lambda x: x[1])
            termNodes.sort(key= lambda x: x[0], reverse=True)
            
            for par in termNodes:
                chi = (par[0]+1, self.treeDepthCnt(par[0]+1)) # depth & index of child
                self.tree[chi] = GlyNode('END', ind=chi[1], parent = self.tree[par])
                    
                    
    
    def drawSbwyMap(self):
        sbwyStps = []
        for k,v in self.tree.items():
            if v.children == ():
                sbwyStps.append(k)
        
        sbwyStps.sort(reverse = True)
        
        for i,stp in enumerate(sbwyStps):
            self.tree[stp].drawSubway(i)

    def buildEncoding(self, maxB, maxC, monoDict, anoDict = {'u': 0, 'a': 1, 'b': 2}):
        """
        

        Parameters
        ----------
        maxB : int
            Max value contained in a subway line in the root nodes of all glycans being considered (maximum number of observed branches).
        maxC : int
            Highest carbon number observed to participate in a glycosidic bond (from all glycans being considered).
        monoDict : dict
            Dictionary linking all observed monosaccharides to a corresponding integer (integers based on monosaccharide frequency rank).
        anoDict : dict, optional
            Dictionary encoding anomeric conformation information as integers.  The default is {'u': 0, 'a': 1, 'b': 2}.

        Returns
        -------
        None, builds encoding in self.encoding as a numpy array, where each row corresponds to a node/monosaccharide and each column corresponds to a descriptor for that node:
            1 -- monosaccharide identity (represented by integer from monoDict)
            2 -- anomeric conformation of the saccharide (0:unknown, 1: alpha, 2: beta)
            3 -- Carbon postions on sacchride participating in glycosidic bonds as a binary converted to a decimal number
                Ex 100001000 (C1 & C6 occupied) --> 264
            4 -- "Subway lines" passing through the node as a binary converted to a decimal number
                Ex 11111100000 (Root node on a glycan with 6 terminal non-reducing saccharides, all 6 subay lines pass through) --> 2016
            5 -- The carbon position of the parent saccharide the node is connected to
            6 -- depth
            7 -- index (differentiate saccharides at the same depth on different branches)

        """
        colNames = ['sugar'] + ['anomeric'] + ['C_links'] + ['B_lines'] + ['parLink', 'sDepth', 'sInd']
        
        self.encoding = np.zeros((len(self.tree.keys()), len(colNames)), dtype = int) # Initialize 2D array to store encoding
        
        for i,nodeKey in enumerate(list(self.tree.keys())):
            base = self.tree[nodeKey] # GlyNode object for the current saccharide
            
            # Prep col 3 value (occupied carbons)
            carbLinks = [str(base.Clinks[str(i)]) for i in range(1,maxC+1)]
            carbLinks = ['1' if c != '0' else '0' for c in carbLinks] # Drop parent/child linked info from each carbon position
            C_binary = ''.join(carbLinks)
            
            # Prep col 4 value (subway lines)
            sbwyLines = ['1' if i in base.subway else '0' for i in range(maxB+1)]
            B_binary = ''.join(sbwyLines)
            
            # Columns 5-7
            liDeIn = [int(base.parentLink), base.depth, base.ind] # link, depth, index info for sugar & parent

            self.encoding[i,] = [monoDict[base.base], anoDict[base.anomeric], binary2decimal(C_binary), binary2decimal(B_binary)] + liDeIn
            
