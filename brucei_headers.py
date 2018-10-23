#!/usr/bin/env python

import sys

def parse(filein):
    prefix = []
    with open(filein, "r+")as input_file:
        print filein
        for line in input_file:
            if line.startswith(">"):
                prefix.append(line.replace(">", "")[0:5])
            else:
                pass
    print set(prefix)




if __name__ == "__main__":
    fileins = sys.argv[1:]
    for i in fileins:
        parse(i)
