from LabelGTCRec import LabelGTC

import os
import sys

import copy

import time

from ..TreeLib import *
from ..TreeLib import TreeUtils, TreeClass


tps1 = time.clock()

s = TreeClass("((A,B),(C,(D,E)));")

print("\n\n\n\n")
print("**********************************************************************************")
print("*                                       INPUT                                    *")
print("**********************************************************************************")
print("\n\n\n\n")

seuil = 0.7

print("------THRESHOLD------")
print("\n")
print(seuil)
print("\n\n\n\n")

print("-----SPECIES TREE-----")
print("\n")
print(s)
print("\n")
print("\n\n\n\n")



g = TreeClass("(((a1_A, b1_B)0, c1_C)0, (((e2_E, e3_E)0, (d2_D, d3_D)0)0, ((d1_D, e1_E)0, c2_C)0)0)0;")

cst = [TreeClass("(a1_A,b1_B);"),TreeClass("c1_C;"), TreeClass("((d1_D, e1_E), c2_C);"), TreeClass("(e2_E, e3_E);"), TreeClass("(d2_D, d3_D);")]

print("-----COVERING SET OF TREES-----")
print("\n")
for tree in cst:
    print(tree)
    print("\n")
print("\n\n\n\n")

lgtc = LabelGTC(s, g, cst, seuil)

print("-----GENES TREE-----")
print("\n")
print(lgtc.getGenesTree().get_ascii(show_internal=True, attributes=["support", "name"]))
print("\n")
print("\n\n\n\n")

print("**********************************************************************************")
print("*                                 ALGORITHM BEGIN                                *")
print("**********************************************************************************")
print("\n\n\n\n")

print("CURRENT TREE BEING ANALIZED")
print(lgtc.getGenesTree().get_ascii(show_internal=True))
print("\n")

lgtc.mergeResolutions()
print("**********************************************************************************")
print("*                            ALGORITHM ENDED SUCCESFULLY                         *")
print("**********************************************************************************")
print("\n\n\n\n")

tps2 = time.clock()
print("Time to compute:")
print(tps2 - tps1)
