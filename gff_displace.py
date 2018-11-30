#!/usr/bin/env python

import sys
from operator import itemgetter

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

## It is important that the gff file should be sorted beforehand.
## Place displacements as (displaced bp, chrom, pos of insertion) in order:
displacements = [(1124, 11, 1767775)]
                #(5110, 14, 768137)]

def displace_gff(sorted_gff_file):
    with open(sorted_gff_file, "r+") as filein:
        for line in filein:
            linelist = line.strip().split("\t")
            if line.startswith("#"):
                #pass
                print line.strip()
            else:
                for d in displacements:
                    if int(chr_ord[linelist[0]]) == int(d[1]):
                        if int(linelist[3]) >= int(d[2]):
                            linelist[3] = str(int(linelist[3]) + d[0])
                            linelist[4] = str(int(linelist[4]) + d[0])

                #pass
                print "\t".join(linelist)


if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		displace_gff(element)
