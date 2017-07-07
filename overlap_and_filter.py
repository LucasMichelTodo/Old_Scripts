#!/usr/bin/env python

import sys
import subprocess as sp
import itertools

if __name__ == "__main__":

	samples = ["10G", "A7K9", "E5K9", "1.2B", "C2"]

	pairs = itertools.permutations(samples, 2)

	for pair in pairs:
		cmd = "peak_overlapper2.py -p1 {} -p2 {}" .format(pair[0]+"_"+pair[1]+"_"+"peaks.xls", pair[0]+"_"+"peaks.xls")
		print cmd
		sp.call(cmd, shell=True)
		







