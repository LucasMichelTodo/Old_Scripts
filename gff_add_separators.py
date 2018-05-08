#!/usr/bin/env python

import sys

def add_separators(gff_file):
    current_chrom = "start"
    with open(gff_file, "r+") as filein:
        for line in filein:
            if line.split()[0] != current_chrom:
                current_chrom = line.split()[0]
                print ""
                print line.strip()
            else:
                print line.strip()




if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		add_separators(element)
