#!/usr/bin/env python
from ..lib.SGT import getMinSGT
from ete3 import Tree


tlist = "((A1__A, B1__B),(C1__C, (C2__C,D2__D)));((C4__C,A2__A),(B2__B,C3__C))"
print (tlist)
gtreelist = [Tree(x+";") for x in tlist.split(';')]
sptree = Tree("((A,B),(C,D));")

#treated_trees = "" #"((C1__C, A1__A),(B1__B, (C2__C,D2__D)));((B1__B, (C2__C,D2__D));(B1__B,C3__C))"

gcontent = "".join([gt.write(format=9) for gt in gtreelist])
scontent = sptree.write(format=9).strip(';')


print "Test with default option ", 
print '------------------------------------------------------'
print getMinSGT(gcontent, scontent, True, "", "")
print "Test with multiple solutions ", 
print '------------------------------------------------------'
print getMinSGT(gcontent, scontent, True, "", "", "", 4)
