#!/usr/bin/env python

import sys
from collections import Counter
import operator
import itertools



def parse_tcd8(fileins):

    all_epitopes = []
    hlas = {}

    for filein in fileins:

        firstline = True
        epitopes = []

        with open(filein, "r+") as infile:
            for line in infile:
                if firstline:
                    #print line.split()[5]
                    firstline = False
                else:
                    linelist = line.split("\t")
                    if float(linelist[7]) <= 0.1:
                        epitopes.append(linelist[5])


        hla = linelist[0]
        ## print hla, set(epitopes)

        hlas[hla] = set(epitopes)

        for i in set(epitopes):
            all_epitopes.append(i)

    counts = Counter(all_epitopes)

    # Sort epitopes by number of counts.
    sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)

    for i in sorted_counts:

        origins = []
        for x in hlas.keys():
            if i[0] in hlas[x]:
                origins.append(x)

        print "{}\t{}\t{}" .format(i[0], i[1], origins)



if __name__ == "__main__":
    filenames = sys.argv[1:]
    parse_tcd8(filenames)
