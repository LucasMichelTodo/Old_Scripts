#!/usr/bin/env python

import sys
import subprocess as sp
from operator import itemgetter

## gt gff3 -sort: Almost does the job bust doesn't sort numerically mRNA/exons/cds from a same gene.

## First we need to fix the gff with: >gt gff3 -tidy file.gff > fixed.gff3
## Then we need to sort it with: >gt gff3 -sort fixed.gff > sorted.gff3
## Then we can proceed with this script

#Transform types to numerical categories to make use of sort.
types_ord = {"gene":1,
            "mRNA":2,
            "exon":3,
            "CDS":4,
            "ncRNA":5,
            "rRNA":5,
            "snoRNA":5,
            "snRNA":5,
            "three_prime_UTR":5,
            "tRNA":5}

def sort_gff(gff_file):
    #cmd = "gt gff3 -tidy {} > {}" .format
    #sp.Popen(cmd, shell=True).wait()
    seqs = []
    with open(gff_file, "r+") as filein:
        header = True
        for line in filein:
            # First skip header lines:
            if header:
                if line.startswith("#"):
                    print line.strip()
                else:
                    header = False
                    seqs.append((line.strip().split())) # First gff line.

            else:
                if line.startswith("###"):
                    for i in seqs:
                        i[2] = types_ord[i[2]]
                    sorted_seqs = sorted(seqs, key=itemgetter(3,2))
                    for i in sorted_seqs:
                        i[2] = [key for key, value in types_ord.iteritems() if value == i[2]][0]
                    for i in sorted_seqs:
                        print "\t".join(i)

                    seqs = []
                else:
                    seqs.append((line.strip().split()))

    for i in seqs:
        i[2] = types_ord[i[2]]
    sorted_seqs = sorted(seqs, key=itemgetter(3,2))
    for i in sorted_seqs:
        i[2] = [key for key, value in types_ord.iteritems() if value == i[2]][0]
    for i in sorted_seqs:
        print "\t".join(i)

if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		sort_gff(element)
