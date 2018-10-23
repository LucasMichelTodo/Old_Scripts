#!/usr/bin/env python

import sys
from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo
from Bio import AlignIO

def get_unmasked_seq(masked_file):
    with open(masked_file, "r+") as infile:
        for line in infile:
            if line.startswith(">"):
                prot = line.replace(">", "").strip()

    clustal_file = masked_file.replace("_consensus00.fasta", "_sorted.aln")
    alignment = AlignIO.read(open(clustal_file), "clustal")

    for record in alignment:
        if record.id == prot:
            print ">"+record.id
            print str(record.seq).replace("-", "")

if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        get_unmasked_seq(i)
