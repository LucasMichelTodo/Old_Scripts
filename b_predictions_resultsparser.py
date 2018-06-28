#!/usr/bin/env python

import sys
import re

seq = ""
seqs = []

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/predictions_antigens.txt", "r+") as infile:
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
            if float(line.split("\t")[2]) >= 0.6:
                seq += line.split("\t")[1]
            else:
                seq += "-"
    #print seq
    seqs.append((seq,antigen))

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/epitopes.fasta", "r+") as infile2:
    for line in infile2:
        if line.startswith(">"):
            pass
        else:
            regex = re.compile(re.escape(str(line.strip())))
            for i in seqs:
                result = re.findall(regex, i[0])
                if result:
                    print result, i[1]
