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

def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def merge_fastas(fasta_list):

    final_fasta = {}
    for element in fasta_list:
        final_fasta.update(fastatodict(element))

    for prot, seq in final_fasta.iteritems():
        if len(prot) > 1:
            print prot.replace(">", ">Exposed_")
            fragments = split_len(seq, 60)
            for fragment in fragments:
                print fragment


if __name__ == "__main__":
    filenames= sys.argv[1:]
    merge_fastas(filenames)
