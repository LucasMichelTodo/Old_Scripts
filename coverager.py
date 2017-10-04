#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
import subprocess
import pysam
import numpy
import csv

sizes = []
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes.append((row[0],row[1]))
		else:
			pass

print sizes

def get_coverage(bamfile):
	bam = pysam.AlignmentFile(bamfile, "rb")
 	coverage = []
 	for chrom in tqdm(sizes):
		cov = bam.count_coverage(reference = chrom[0], start = 1, end = int(chrom[1]))
		total_cov = numpy.sum(cov, axis=0)
		coverage.append(total_cov) ## Assuming coding regions do NOT OVERLAP
	
	coverage = numpy.concatenate(coverage, axis=0)
	print coverage
	with open(element.replace(".bam", "_coverage.txt"), "a") as results_file:
		results_file.write(coverage)

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		get_coverage(element)


