#!/usr/bin/env python

import sys
from collections import Counter
import operator



def cross_epitopes(binders_list):

    with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_27_02_18/fastas/tcruzi_all_cons_rankpep_pedro.txt", "r+") as filein1:
        pedro_epitopes = []
        for line in filein1:
            pedro_epitopes.append(line.split()[0])

    for filein in binders_list:
        with open(filein, "r+") as infile:
            epitopes = []
            for line in infile:
                epitopes.append(line.split()[2])

    set1 = set(pedro_epitopes)
    set2 = set(epitopes)
    print len(set1)
    print len(set2)
    print len(set2.intersection(set1))

if __name__ == "__main__":
    cross_epitopes(sys.argv[1:])
