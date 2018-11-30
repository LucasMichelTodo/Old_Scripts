#!/usr/bin/env python

ags = []
with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_proteomes.fasta") as infile:
    for line in infile:
        if line.startswith(">"):
            ags.append(line.strip()[0:3])

print set(ags)
print len(set(ags))
