#!/usr/bin/env python

import sys
from Bio import SeqIO


# Convert a big fasta file into many smaller fastas with nseq sequences each:
def divide_fasta(fasta_file, nseq):
    current = []
    counter = 0
    filenum = 1
    for record in SeqIO.parse(fasta_file, "fasta"):

        if counter < nseq:
            counter += 1
            current.append(record)

            if counter == nseq:
                SeqIO.write(current, "fasta_file_{}.fasta" .format(filenum), "fasta")
                counter = 0
                current = []
                filenum += 1


if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        divide_fasta(file, 1)
