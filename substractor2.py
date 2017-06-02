#!/usr/bin/env python

import pysam
import csv 
import sys
import numpy
import pandas as pd



depths = {"10G":4535787, "1.2B":4997472, "A7": 10915415, "C2": 5176848, "E5": 9366386}

def coverage(bamfile):

	sam = pysam.AlignmentFile(bamfile, "rb")

	chrom_cov = {}
	sizes = {}
	
	with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
		chrom_sizes = csv.reader(csvfile, delimiter='\t')
		for row in chrom_sizes:
			if len(row) == 2:
				sizes[row[0]] = row[1]
			else:
				pass

	for key in sizes:
		coverage = []
		cov = sam.count_coverage(reference=key, start=1, end=int(sizes[key]))
		total_cov = numpy.sum(cov, axis=0)
		smooth_cov = pd.Series(total_cov).rolling(2000).mean().tolist()
		chrom_cov[key] = smooth_cov

	return chrom_cov

def normalize_coverage(bamfile):
	
	norm_cov = {}
	abs_cov = coverage(bamfile)
	for key in depths:
		if key in bamfile:
			for chrom in abs_cov:
				norm_cov[chrom] = abs_cov[chrom]/float(depths[key])*1000000
				 
	return norm_cov




def susbstract_norm_cov(bam1, bam2):

	sub_cov = {}
	dict1 = normalize_coverage(bam1)
	dict2 = normalize_coverage(bam2)

	for key in dict1:
		sub_cov[key] = dict1[key] - dict2[key]

	for key in sub_cov:
		i = 0
		for element in sub_cov[key]:
			with open("result_2.txt", "a+") as bed_file:
				bed_file.write("{}\t{}\t{}\t{}\n" .format(key, i, i+1, element))
				i += 1

	return sub_cov




if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	susbstract_norm_cov(filenames[0], filenames[1])