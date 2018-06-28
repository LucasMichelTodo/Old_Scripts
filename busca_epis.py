#!/usr/bin/env python

i = 1
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/New_strategy/cluster_withref/msa/sorted_epitopes.txt", "r+") as filein:
    for line in filein:
        exec "prot = "+line
        if prot[1] == 0.0:
            print ">ep"+str(i)+"_"+prot[2]
            print prot[0]
            i += 1
