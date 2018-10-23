#!/usr/bin/env python

import sys

epitopes = []
exposed_seqs = []
seq = ""

with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/data/tritrypdb_output/aat_vsg_2_fastas.txt", "r+") as infile:
#with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/data/tritrypdb_output/aat_isg_2_fastas.txt", "r+") as infile:
    for line in infile:
        if line.startswith(">"):
            exposed_seqs.append(seq)
            seq = ""
        else:
            seq += line.strip()

exposed_seqs.append(seq)


with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_consensus00_onlyAA_15aa.fasta", "r+") as infile1:
    for line in infile1:
        if line.startswith(">"):
            pass
        else:
            epitopes.append(line.strip())

#print epitopes
#print exposed_seqs

for epitope in epitopes:
    for seq in exposed_seqs:
        if epitope in seq:
            print epitope
            print seq
            print "------------------------"
