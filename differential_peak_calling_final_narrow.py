#!/usr/bin/env python

import sys
import subprocess as sp
import itertools

if __name__ == "__main__":

	samples = ["10G", "A7K9", "E5K9", "1.2B", "C2"]
	for sample in samples:
		cmd = "macs2 callpeak -g 2.3e+7 -f BAMPE -t {}_me_sort_q5.bam -c {}_in_sort_q5.bam -n {} --nomodel --extsize 117 --fe-cutoff 3" .format(sample, sample, sample)
		print cmd
		sp.call(cmd, shell=True)

	pairs = itertools.permutations(samples, 2)

	for pair in pairs:
		cmd = "macs2 callpeak -g 2.3e+7 -f BAMPE -t {}_me_sort_q5.bam -c {}_me_sort_q5.bam -n {} --nomodel --extsize 117 --fe-cutoff 3" .format(pair[0], pair[1], pair[0]+"_"+pair[1])
		print cmd
		sp.call(cmd, shell=True)

	for pair in pairs:
		cmd = "peak_overlapper2.py -p1 {} -p2 {}" .format(pair[0]+"_"+pair[1]+"_"+"peaks.xls", pair[0]+"_"+"peaks.xls")
		print cmd
		sp.call(cmd, shell=True)
		







