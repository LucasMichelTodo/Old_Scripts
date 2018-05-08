#!/usr/bin/env python

import sys

def to_bed(table_file):
    with open(table_file, "r+") as file1:
        header = True
        for line in file1:
            if header:
                header = False
            else:
                line_list = line.strip().split("\t")
                line_list.insert(2, line_list[1])
                print "\t".join(line_list)


if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
	       to_bed(element)
