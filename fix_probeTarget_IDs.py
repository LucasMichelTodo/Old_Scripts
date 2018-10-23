#!/usr/bin/env python

import sys

def translate(filein):
	with open(filein) as file:
		for line in file:
			genes.append(line.strip().split(":")[0])



if __name__ == "__main__":
	filein = sys.argv[1]
	translate(filein)
