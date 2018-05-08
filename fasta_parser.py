#!/usr/bin/env python

import sys

def merge_fastas(fasta_file):
    with open(fasta_file, "r+") as file1:
        for line in file1:
            if line.startswith(">"):
                print line.split(" | ")


if __name__ == "__main__":
    filenames= sys.argv[1:]
    print filenames
    for filein in filenames:
        merge_fastas(filein)
