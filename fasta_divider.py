#!/usr/bin/env python

import sys
from Bio import SeqIO


# Convert a big fasta file into many smaller fastas with nseq sequences each:
def divide_fasta(fasta_file, nseq):
    current_100 = []
    counter = 0
    filenum = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        if counter < 99:
            counter += 1
            current_100.append(record)
        else:
            SeqIO.write(current_100, "hundred_{}.fasta" .format(filenum), "fasta")
            counter = 0
            current_100 = []
            filenum += 1

if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        divide_fasta(file, 100)
