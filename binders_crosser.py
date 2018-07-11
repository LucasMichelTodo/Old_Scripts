#!/usr/bin/env python

import sys
from collections import Counter
import operator
import itertools


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

    hla_rank = {
    "HLA_A0201":21.23,
    "HLA_A0101":10.16,
    "HLA_A0202":0.41,
    "HLA_A0203":0.0,
    "HLA_A0205":1.73,
    "HLA_A0206":0.0,
    "HLA_A0207":0.0,
    "HLA_A6801":2.48,
    "HLA_A6802":1.39,
    "HLA_A0301":7.57,
    "HLA_A1101":5.71,
    "HLA_A3101":6.01,
    "HLA_A3301":1.87,
    "HLA_A6601":0.61,
    "HLA_A2402":9.58,
    "HLA_B3801":3.7,
    "HLA_B0702":5.16,
    "HLA_B3501":6.18,
    "HLA_B5101":7.88,
    "HLA_B5102":0.03,
    "HLA_B5301":0.95,
    "HLA_B5401":0.0,
    "HLA_B1501":2.41,
    "HLA_B1502":0.0,
    }

    epitope_HLAs = {}

    # Iterate over files and populate a dictionary with epitopes as keys and a list of HLAs as value:
    for filein in binders_list:
        with open(filein, "r+") as infile:
            for line in infile:
                epitope_HLAs.setdefault(line.split()[2], []).append(filein.replace("./", "").replace("_binders_SB_filetered.txt", ""))

    # Define a function to calculate coverage from list of HLAs:
    def calculate_rank(hla_list):
        coverage = 0
        for hla in set(hla_list):
            coverage += hla_rank[hla]
        return coverage

    #Using the created function, sort epitopes by coverage:
    sorted_counts = sorted(epitope_HLAs.items(), key=lambda x: calculate_rank(x[1]), reverse=True)

    # Print results
    # Header
    print "Epitope\tCoverage\tn_HLAs\tn_total\tHLAs"

    for i in range(0,3):
        a = i
        print sorted_counts[i]
        target = 95
        current_epilsit = [sorted_counts[i][0]]
        current_hlalist = sorted_counts[i][1]
        current_sum = calculate_rank(current_hlalist)
        while current_sum < target:
            if calculate_rank(current_hlalist+sorted_counts[a][1]) > current_sum:
                #print current_hlalist
                #print sorted_counts[a][1]
                current_epilsit.append(sorted_counts[a][0])
                current_hlalist.extend(sorted_counts[a][1])
                current_sum = calculate_rank(current_hlalist)
                print current_sum
                a += 1
            else:
                a += 1
        print current_epilsit
        print set(current_hlalist)


        set_test = ["poma", "pera", "platan"]
        combs = itertools.combinations(set_test, 2)
        for i in combs:
            print i
        print combs

        # else:
        #     print current_epilsit
        #     print current_sum


        # while current_sum < target:
        #     for x in sorted_counts:
        #         current_epilsit.append(x[0])
        #         current_hlalist.extend(x[1])
        #         current_sum = calculate_rank(current_hlalist)
        #         if current_sum >= target:
        #             print current_epilsit
        #             print current_sum


    # for i in sorted_counts:
    #     # Count total times the epitope appears:
    #     hla_counts = Counter(i[1])
    #
    #     # Format string with each HLA counts sepparately:
    #     hla_counts_string = ""
    #     for hla, count in hla_counts.items():
    #         hla_counts_string += "{}({})," .format (hla, count)
    #     hla_counts_string = hla_counts_string.strip(",")
    #
    #     # Result body:
    #     print "{}\t{}\t{}\t{}\t{}" .format(i[0], calculate_rank(i[1]), len(set(i[1])), len(i[1]), hla_counts_string)


if __name__ == "__main__":
    cross_epitopes(sys.argv[1:])
