#!/usr/bin/env python

import sys
from Bio import SeqIO
input_file = open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/proteomes_orfs_aa/all_fasta1.fasta")

ref_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

# for key, value in ref_dict.iteritems():
#     print key
#     print value


def recover_fasta(ag_file):
    prots = []
    with open(ag_file, "r+") as infile:
        for line in infile:
            prots.append(line.strip())


    for i in set(prots):
        print ">{}" .format(i)
        try:
            print ref_dict[i].seq
        except:
            print ref_dict[i+":mRNA-p1"].seq


if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		recover_fasta(file)
