"""Testing Multi-PolyRes use"""

from LabelGTCRec import LabelGTC

import os
import sys

import copy

import time

from ..TreeLib import *
from ..TreeLib import TreeUtils, TreeClass


tps1 = time.clock()

s = TreeClass("((A,B),C);")

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



g = TreeClass("((((a_A,x_B)0.2,(b_B,e_C)0.2)0.2,y_C)0.2,((i_B,k_A)0.1,((c_C, j_A)0.1,(d_C,(g_A,h_A)0.2)0.8)0.2)0.8)0.2;")

cst = [TreeClass("a_A;"),TreeClass("x_B;"), TreeClass("y_C;"), TreeClass("b_B;") , TreeClass("e_C;"), TreeClass("c_C;") ,TreeClass("j_A;"), TreeClass("g_A;") ,TreeClass("h_A;"), TreeClass("d_C;"), TreeClass("i_B;"), TreeClass("k_A;")]

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
