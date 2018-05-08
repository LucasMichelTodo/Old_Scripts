#!/usr/bin/env python

import sys
from collections import Counter
import operator


## Cross and keep only elements appearing in every list: ###############################################################################################

# def cross_epitopes(binders_list):
#     haplotypes = []
#     for filein in binders_list:
#         with open(filein, "r+") as infile:
#             epitopes = []
#             for line in infile:
#                 epitopes.append(line.split()[2])
#             haplotypes.append(epitopes)
#
#     result = set(haplotypes[0])
#     for s in haplotypes[1:]:
#         result.intersection_update(s)
#         print len(result)
#
#     for i in result:
#         print i

## Cross lists and create a counter for each epitope with number of appearences. #######################################################################

# def cross_epitopes(binders_list):
#     # Create empty list for epitopes:
#     epitopes = []
#
#     # Open every file and dump epitopes in the list.
#     for filein in binders_list:
#         with open(filein, "r+") as infile:
#             for line in infile:
#                 epitopes.append(line.split()[2])
#
#     # For each epitope, count number of appearences in the list.
#     counts = Counter(epitopes)
#
#     # Sort epitopes by number of counts.
#     sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
#
#     for i in sorted_counts:
#         print "{}   {}" .format(i[0], i[1])

## Create file with epitope, number of appearences and HLAs ############################################################################################

def cross_epitopes(binders_list):

    epitope_HLAs = {}

    # Iterate over files and populate a dictionary with epitopes as keys and a set of HLAs as value:
    for filein in binders_list:
        with open(filein, "r+") as infile:
            for line in infile:
                epitope_HLAs.setdefault(line.split()[2], []).append(filein.replace("./", "").replace("_binders_SB_filetered.txt", ""))

    # Create a list of sorted tuples from the dictionary. Sort by lengths of the HLAs set.
    sorted_counts = sorted(epitope_HLAs.items(), key=lambda x: len(set(x[1])), reverse=True)

    # Print results
    # Header
    print "Epitope\tn_HLAs\tn_total\tHLAs"

    for i in sorted_counts:
        # Count total times the epitope appears:
        hla_counts = Counter(i[1])

        # Format string with each HLA counts sepparately:
        hla_counts_string = ""
        for hla, count in hla_counts.items():
            hla_counts_string += "{}({})," .format (hla, count)
        hla_counts_string = hla_counts_string.strip(",")

        # Result body:
        print "{}\t{}\t{}\t{}" .format(i[0], len(set(i[1])), len(i[1]), hla_counts_string)


if __name__ == "__main__":
    cross_epitopes(sys.argv[1:])
