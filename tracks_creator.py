#!/usr/bin/env python

import pysam
import csv
import sys
import numpy
import pandas as pd


# #Numberof reads after filtering in treatment:
# depths = {"10G":4535787.0, "1.2B":4997472.0, "A7": 10915415.0, "C2": 5176848.0, "E5": 9366386.0}

#Create a dictinary with size of each chromosome and a lsit of ordered chromosomes:
sizes = {}
chroms=[]

with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
			chroms.append(row[0])
		else:
			pass

#Calculate coverage over whole genome:
def coverage(bamfile):

	sam = pysam.AlignmentFile(bamfile, "rb")

	chrom_cov = {}

	for key in sizes:
		coverage = []
		cov = sam.count_coverage(reference=key, start=1, end=int(sizes[key]))
		total_cov = numpy.sum(cov, axis=0)
		#print "total_cov"
		#print total_cov[1:10]
		chrom_cov[key] = total_cov

	return chrom_cov

# "Normalize" coverage diving by reads after filtering and multiplying by 1000000:
def rpm_coverage(bamfile):

	# Get number of mapped reads of the alignment.
	sam = pysam.AlignmentFile(bamfile, "rb")
	depth = sam.mapped

	norm_cov = {}
	abs_cov = coverage(bamfile)
	for chrom in abs_cov:
		norm_cov[chrom] = abs_cov[chrom]/float(depth)*1000000
		#print "norm_cov"
		#print norm_cov[chrom][1:10]

	return norm_cov


# Substract two "normalized" coverage traks:
def susbstract_rpm_cov(bam1, bam2):

	sub_cov = {}
	dict1 = rpm_coverage(bam1)
	dict2 = rpm_coverage(bam2)

	for key in dict1:
		sub_cov[key] = dict1[key] - dict2[key]

	for chrom in chroms:
		i = 0
		for base_cov in sub_cov[chrom]:
			print "{}\t{}\t{}\t{}" .format(chrom, i, i+1, base_cov)
			i += 1

	return sub_cov


def ratio_norm_cov(bam1, bam2):

	ratio_cov = {}
	dict1 = normalize_coverage(bam1)
	dict2 = normalize_coverage(bam2)

	for key in dict1:
		ratio_cov[key] = numpy.log(dict1[key] / dict2[key])

	for chrom in chroms:
		i = 0
		for base_cov in ratio_cov[chrom]:
			print "{}\t{}\t{}\t{}" .format(chrom, i, i+1, base_cov)
			i += 1


	return ratio_cov



def rpm_normby_input(treat, control):

	input_cov = {}
	treat_cov = rpm_coverage(treat)
	treat_cov_pseudo = {}
	for key in treat_cov:
		treat_cov_pseudo[key] = treat_cov[key]+1

	control_cov = rpm_coverage(control)
	control_cov_pseudo = {}
	for key in control_cov:
		control_cov_pseudo[key] = control_cov[key]+1

	for key in treat_cov_pseudo:
		input_cov[key] = numpy.log2((treat_cov_pseudo[key]/control_cov_pseudo[key]))

	for chrom in chroms:
		i = 0
		for base_cov in input_cov[chrom]:
			print "{}\t{}\t{}\t{}" .format(chrom, i, i+1, base_cov)
			i += 1


	return input_cov


def print_coverage_at_intervals(bamfile, inter):

	dict1 = normalize_coverage(bamfile)

	for key in dict1:
		pos = 0
		i = 0
		cov = 0

		for element in dict1[key]:
			cov += element
			i += 1
			pos += 1

			if i == int(inter):

				with open(bamfile.replace(".bam", "_fullcoverage.bed"), "a+") as bed_file:
					print bamfile.replace(".bam", "_fullcoverage.bed")
					bed_file.write("{}\t{}\t{}\t{}\n" .format(key, pos-i, pos, cov))
					i = 0
					cov = 0
					pos = pos+i

###### Change this part to use whatever function you like.

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	rpm_normby_input(filenames[0], filenames[1])
