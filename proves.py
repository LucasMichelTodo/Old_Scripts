#!/usr/bin/env python

import sys
import pysam
import numpy
import csv
import pandas as pd

depths = {"10G":4535787, "1.2B":4997472, "A7K9": 10915415, "C2": 5176848, "E5K9": 9366386} # comprovar que corespon a q5

def get_coverage(bamfile, chrom, start, stop):
	bam = pysam.AlignmentFile(bamfile, "rb")
	cov = bam.count_coverage(reference=chrom, start=start, end=stop)
	print cov
	total_cov_vect = numpy.sum(cov, axis=0)
	print total_cov_vect
	total_cov = sum(total_cov_vect) / float(stop-start)
	print total_cov

	return total_cov

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		get_coverage(element)