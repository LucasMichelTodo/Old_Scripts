#!/usr/bin/env python

import sys
from operator import itemgetter

#Transform types to numerical categories to make use of sort.
types_ord = {"gene":1,
            "mRNA":2,
            "exon":8,
            "CDS":9,
            "ncRNA":3,
            "rRNA":4,
            "snoRNA":5,
            "snRNA":6,
            "three_prime_UTR":10,
            "tRNA":7}

chr_ord = {"Pf3D7_01_v3":1,
        "Pf3D7_02_v3":2,
        "Pf3D7_03_v3":3,
        "Pf3D7_04_v3":4,
        "Pf3D7_05_v3":5,
        "Pf3D7_06_v3":6,
        "Pf3D7_07_v3":7,
        "Pf3D7_08_v3":8,
        "Pf3D7_09_v3":9,
        "Pf3D7_10_v3":10,
        "Pf3D7_11_v3":11,
        "Pf3D7_12_v3":12,
        "Pf3D7_13_v3":13,
        "Pf3D7_14_v3":14,
        "Pf3D7_API_v3":15,
        "Pf_M76611":16}

def sort_gff(gff_file):
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
                seqs.append((line.strip().split()))


    for i in seqs:
        i[0] = chr_ord[i[0]]
        i[2] = types_ord[i[2]]
        i[3] = int(i[3])

    sorted_seqs = sorted(seqs, key=itemgetter(0,3,2))
    for i in sorted_seqs:
        i[0] = [key for key, value in chr_ord.iteritems() if value == i[0]][0]
        i[2] = [key for key, value in types_ord.iteritems() if value == i[2]][0]
        i[3] = str(i[3])
    for i in sorted_seqs:
        print "\t".join(i)

if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		sort_gff(element)
