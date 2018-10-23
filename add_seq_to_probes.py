#!/usr/bin/env python

# Import packages
import sys

probes = []

def add_seq_to_probelist(probe_list):
    with open(probe_list, "r+") as infile1:
        for line in infile1:
            probes.append(line.strip())

    d = {}
    with open("/home/lucas/ISGlobal/Arrays/Array_Annotation/probes.fasta", "r+") as infile2:
        for line in infile2:
            if line.startswith(">"):
                current_probe = line.strip().replace(">","").replace(":", "/")
            else:
                d[current_probe] = line.strip()

    for i in probes:
        print i+"\t"+d[i]

if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		add_seq_to_probelist(element)
