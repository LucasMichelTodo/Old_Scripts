#!/usr/bin/env python

import sys
import re

def calculate_stats(fasta_file):

    aa = re.compile("[A-Z]")
    gap = re.compile("-")
    masked = re.compile("\*")

    aa_counter = 0
    gap_counter = 0
    masked_counter = 0

    with open(fasta_file, "r+") as filein:
        for line in filein:
            if line.startswith(">"):
                pass
            else:
                aa_counter += len(re.findall(aa, line))
                gap_counter += len(re.findall(gap, line))
                masked_counter += len(re.findall(masked, line))

    print "Number of conserved residues {} ({}%)" .format(aa_counter, round(float(aa_counter)*100/(aa_counter+masked_counter+gap_counter),2))
    print "Number of masked residues {} ({}%)" .format(masked_counter, round(float(masked_counter)*100/(aa_counter+masked_counter+gap_counter),2))
    print "Number of gaps {} ({}%)" .format(gap_counter, round(float(gap_counter)*100/(aa_counter+masked_counter+gap_counter),2))

if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        calculate_stats(file)
