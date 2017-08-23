#!/usr/bin/env python
from minSGT import getMinSGT
from ete3 import Tree


tlist = "((A1__A, C1__C),(B1__B, (C2__C,D2__D)));((B2__B,D3__D),(B3__B,C3__C))"
print (tlist)
gtreelist = [Tree(x+";") for x in tlist.split(';')]
sptree = Tree("((A,B),(C,D));")

#treated_trees = "" #"((C1__C, A1__A),(B1__B, (C2__C,D2__D)));((B1__B, (C2__C,D2__D));(B1__B,C3__C))"

clades_to_preserve = "(C2__C:1,D2__D:1)0.1:1;"
treated_trees = "((C1__C, A1__A),(B1__B, (C2__C,D2__D)));((B1__B, (C2__C,D2__D));(B1__B,C3__C))"

gcontent = "".join([gt.write(format=9) for gt in gtreelist])
scontent = sptree.write(format=9).strip(';')

print "Test with clade to preserve: ", clades_to_preserve
print '------------------------------------------------------'
print getMinSGT(gcontent, scontent, True, clades_to_preserve, treated_trees, "")
print "\nTest without clade to preserve"
print '------------------------------------------------------'

print getMinSGT(gcontent, scontent, True, "", treated_trees, "")
