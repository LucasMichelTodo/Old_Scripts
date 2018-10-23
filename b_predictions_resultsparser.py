#!/usr/bin/env python

import sys
import re

seq = ""
seqs = []

def parse_bepipred2(filein):

    with open(filein, "r+") as infile:
        for line in infile:
            if line.startswith("input:"):
                #print seq
                try:
                    seqs.append((seq,antigen))
                except:
                    pass
                seq = ""
                antigen = line.strip().split()[1]
                #print ">{}" .format(antigen)
            elif line.startswith("Position"):
                pass
            else:
                if float(line.split("\t")[2]) >= 0.5:
                    seq += line.split("\t")[1]
                else:
                    seq += "-"
        #print seq
        seqs.append((seq,antigen))


    count = 0
    for i in seqs:
        print ">"+i[1]
        print i[0]

    
if __name__ == "__main__":
    filename= sys.argv[1]
    parse_bepipred2(filename)
