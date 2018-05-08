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

def merge_fastas(fasta_list):

    final_fasta = {}
    for element in fasta_list:
        final_fasta.update(fastatodict(element))

    for prot, seq in final_fasta.iteritems():
        if len(prot) > 1:
            print ">exposed | "+prot
            print seq


if __name__ == "__main__":
    filenames= sys.argv[1:]
    merge_fastas(filenames)
