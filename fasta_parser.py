#!/usr/bin/env python

tags = []
esme = []

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/proteomes_orfs_aa/all_fasta_noMar.fasta", "r+") as file1:
    for line in file1:
        if line.startswith(">"):
            tags.append(line[0:6])

        if line.startswith(">"):
            split = line.split(" | ")
            for element in split:
                if element.startswith("organism"):
                    esme.append((element, split[0][0:6]))


for i in set(tags):
    print i

for i in set(esme):
    print i
