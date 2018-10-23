#!/usr/bin/env python

import sys

def remove_introns(genomic_seq):

    seq = ""
    with open(genomic_seq, "r+") as infile:
        for line in infile:
            if line.startswith(">"):
                ag = ">lncasGDV1 | "+line.strip().replace(">", "")
            else:
                seq += line.strip()

    print ag
    print seq[0:405]+seq[687:727]+seq[926:2728]+seq[3574:3824]+seq[4316:4739]



if __name__ == "__main__":
    filename= sys.argv[1]
    remove_introns(filename)
