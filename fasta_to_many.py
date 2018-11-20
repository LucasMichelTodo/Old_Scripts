#!/usr/bin/env python

import sys

def fastatomany(fasta_file):

    ## First convert the input fasta into a dict to iterate over it.
    fasta = {}
    seq = ""
    prot = ""

    with open(fasta_file, "r+") as file1:
        for line in file1:
            if line.startswith(">"):
                fasta[prot] = seq
                seq = ""
                prot = line.strip()
            else:
                seq += line.strip()

                fasta[prot] = seq


    ## For item in the dict write a fasta file
    i = 1
    for record, seq in fasta.iteritems():
        #print i
        with open("fasta_{}.fasta" .format(i), "w+") as current_fasta:
            current_fasta.write(record+"\n")
            current_fasta.write(seq)
            i+=1



if __name__ == "__main__":
    filenames= sys.argv[1]
    fastatomany(filenames)
