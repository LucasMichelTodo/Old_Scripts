#!/usr/bin/env python

import sys


def fastatodict(fasta_file):

    fasta = {}
    seq = ""
    prot = ""

    with open(fasta_file, "r+") as file1:
        for line in file1:
            if line.startswith(">"):
                fasta[prot] = seq
                seq = ""
                prot = line.strip()
            else:
                seq += line.strip()

                fasta[prot] = seq

    return fasta


ref_fasta = fastatodict("/home/lucas/ISGlobal/Chip_Seq_Oriol/Chr14/chr14.fasta") 

mod_fasta = fastatodict("/home/lucas/ISGlobal/Chip_Seq_Oriol/Reference_Genomes/PlasmoDB-39_Pfalciparum3D7_Genome_Ind_System_In_Revised.fas")


print "Reference Fasta:"
for key, value in ref_fasta.iteritems():
    print key, value[-10:-1]

print "Modified Fasta:"
for key, value in mod_fasta.iteritems():
    print key, len(value)


# print ">Reference_Chr14"
# print ref_fasta[">Pf3D7_14_v3"][762614:773407]
# print ">Inducible_Chr14"
# print mod_fasta[">Pf3D7_14_v3"][762614:773407]

