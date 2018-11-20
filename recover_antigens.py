#!/usr/bin/env python

antigens = []
with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_consensus00_onlyAA_15aa.fasta", "r+") as file1:
    for line in file1:
      if line.startswith(">"):
          prot = line.split("/")[0].strip()
          antigens.append(prot)

fasta = {}
seq = ""
prot = ""

with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_proteomes.fasta", "r+") as file2:
    for line in file2:
        if line.startswith(">"):
            fasta[prot] = seq
            seq = ""
            prot = line.split(" | ")[0].strip()
        else:
            seq += line.strip()

    fasta[prot] = seq

for ag in antigens:
    if ag in fasta.keys():
        print ag
        print fasta[ag]
        del fasta[ag]
