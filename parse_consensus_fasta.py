#!/usr/bin/env python

import sys

# This script removes gaps ("-") and masked AA ("*") and makes a new entry for each resulting fragment in a new fasta file.
def parse_consensus(consensus_file):

    # Convert the fasta file into a dict.
    with open(consensus_file, "r+") as file1:
        fasta = {}
        for line in file1:
            if line.startswith(">"):
                prot = line.strip()
                fasta[prot] = ""
            else:
                fasta[prot] += line.strip()
        #print fasta

    for prot, seq in fasta.iteritems():
        fragment = ""
        fragments = []
        for i in seq:
            if i not in "-*":
                fragment += i
            else:
                if len(fragment) > 8:
                    fragments.append(fragment)
                fragment = ""
        if len(fragment) > 8:
            fragments.append(fragment)


        if len(fragments) >= 1:
            i = 1
            for x in fragments:
                print "{}_{}" .format(prot, i)
                print x
                i += 1



if __name__ == "__main__":
    filenames= sys.argv[1:]
    for filein in filenames:
        parse_consensus(filein)
