#!/usr/bin/env python

import sys
from collections import Counter
import operator
import itertools



def parse_tcd4(fileins):

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
                    linelist = line.split()
                    if float(linelist[6]) <= 0.01:
                        if linelist[5] == "COMB.LIB.-SMM-NN":
                            epitopes.append(linelist[7])
                        elif linelist[5] == "SMM-NN-Sturniolo":
                            epitopes.append(linelist[10])
                        elif linelist[5] == "NetMHCIIpan":
                            epitopes.append(linelist[16])
                        elif linelist[5] == "SMM":
                            epitopes.append(linelist[10])
                        else:
                            print "Et falta!", filein, linelist[5]

        hla = filein.replace("_predictions.txt", "")
        #print hla, set(epitopes)
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
    parse_tcd4(filenames)
