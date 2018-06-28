#!/usr/bin/env python

import sys
from Bio import SeqIO
# input_file = open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/proteomes_orfs_aa/all_fasta1.fasta")
#
# ref_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

# for key, value in ref_dict.iteritems():
#     print key
#     print value

prots = []
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/epitopes.fasta", "r+") as infile:
    for line in infile:
        if line.startswith(">"):
            if line.split("_")[1] == "Exposed":
                #prots.append(line.split("_")[2].strip().split("/")[0])
                prot = line.split("_")[2].strip().split("/")[0]
            else:
                prot = line.split("_")[1].strip().split("/")[0]
                #prots.append(line.split("_")[1].strip().split("/")[0])
        else:
            epitope = line.strip()
            print epitope, prot


# for i in set(prots):
#     print ">{}" .format(i)
#     try:
#         print ref_dict[i].seq
#     except:
#         print ref_dict[i+":mRNA-p1"].seq
