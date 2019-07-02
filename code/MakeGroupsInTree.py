#!/usr/bin/python

import os
import numpy as np
from ete3 import Tree
import sys

treefile = sys.argv[1]
N = int(sys.argv[2])
## STEP 1: read tree
T = Tree(treefile)
node2leaves = T.get_cached_content()
T.NbOfGroups = N
GROUPS = []
## STEP 2: traverse the tree and store the number of descendants of each node
for n in T.traverse(): 
	n.nbdesc = len(node2leaves[n])
	for u in n.get_children():
		if (n.NbOfGroups>1):
			PropForThisChild = len(node2leaves[u])/len(node2leaves[n])
			u.NbOfGroups = round(PropForThisChild*n.NbOfGroups)
		else: 
			u.NbOfGroups = 0
	if (n.NbOfGroups==1):
		GROUPS.append(n.get_leaf_names())

cpt = 0
for k in GROUPS:
	cpt = cpt+1
	for sp in k:
		print(sp,"\t",cpt)

