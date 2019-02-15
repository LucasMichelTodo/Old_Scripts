#!/usr/bin/env python

import sys


def insert_tag(fasta_file, tag):
    with open(fasta_file, "r+") as infile:
        for line in infile:
            if line.startswith(">"):
                print line.replace(">", ">"+tag+"_").strip()
            else:
                print line


if __name__ == "__main__":
    filenames = sys.argv[1:]
    insert_tag(filenames[0], filenames[1]).strip()
