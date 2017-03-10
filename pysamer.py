#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
from Bio.Seq import Seq
import subprocess
import pysam
from tqdm import tqdm
import csv
import pandas
import numpy
from Bio.Seq import Seq
from pf3D7_coding_regions import coding, noncoding


### FUNCTIONS

# Get fragment lengths for all reads and print to a csv file.

def get_frag_len(samfile):
	if samfile[-4:] == ".sam":
		sam = pysam.AlignmentFile(samfile, "r")
	elif samfile[-4:] == ".bam":
		sam = pysam.AlignmentFile(samfile, "rb")
	else:
		print "Not a bam nor a sam file!"
		

	lens = []
	for read in sam.fetch():
		lens.append(read.template_length)
	with open(samfile.replace(samfile[-4:], "_lengths.csv"),"w") as csvfile:
		writer = csv.writer(csvfile, delimiter = "\t", quoting=csv.QUOTE_MINIMAL)
		writer.writerow(lens)
	sam.close()

# Get MAPQ for all reads and print to a csv file:

def get_MAPQ(samfile):
	if samfile[-4:] == ".sam":
		sam = pysam.AlignmentFile(samfile, "r")
	elif samfile[-4:] == ".bam":
		sam = pysam.AlignmentFile(samfile, "rb")
	else:
		print "Not a bam nor a sam file!"
		

	mapq = []
	for read in sam.fetch():
		mapq.append(read.mapping_quality)
	with open(samfile.replace(samfile[-4:], "_MAPQ.csv"), "w") as csvfile:
		writer = csv.writer(csvfile, delimiter = "\t", quoting=csv.QUOTE_MINIMAL)
		writer.writerow(mapq)
	sam.close()

#Get unaligned reads via samtools:

def get_unaligned(samfile):
	cmd = "samtools view -h -f 4 {} > {}" .format(samfile, samfile.replace(".sam", "_unmapped.sam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools view -h -F4 -f8 {} > {}" .format(samfile, samfile.replace(".sam", "_mateunmapped.sam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools merge -h {} {} {} {}" .format(samfile, samfile.replace(".sam","_somethingununmapped.sam"), 
												samfile.replace(".sam", "_unmapped.sam"),
												samfile.replace(".sam", "_mateunmapped.sam"))
	subprocess.call(cmd, shell=True)

#Get number of reads with an unaligned mate:
	
def get_unaligned_mate(samfile):
	if samfile[-4:] == ".sam":
		sam = pysam.AlignmentFile(samfile, "r")
	elif samfile[-4:] == ".bam":
		sam = pysam.AlignmentFile(samfile, "rb")
	else:
		print "Not a bam nor a sam file!"
		
	unpaired_mates = 0
	for read in sam.fetch():
		if read.mate_is_unmapped:
			unpaired_mates += 1
	print "Unpaired mates: {}" .format(unpaired_mates)
	sam.close()

#Print flags prensent in alignment (type and amount):

def get_flags(samfile):
	if samfile[-4:] == ".sam":
		sam = pysam.AlignmentFile(samfile, "r")
	elif samfile[-4:] == ".bam":
		sam = pysam.AlignmentFile(samfile, "rb")
	else:
		print "Not a bam nor a sam file!"
		
	flags = {}
	for read in sam.fetch():
		if read.flag in flags:
			flags[read.flag] += 1
		else:
			flags[read.flag] = 1
	print "Flags present: {}" .format(flags)
	sam.close()

#Get coverage splitted among coding and noncoding regions, written into a .txt file:

def get_coding_coverage(bamfile):
	sam = pysam.AlignmentFile(bamfile, "rb")
	
	references = []
	starts = []
	ends = []
	for key in coding:
		for i in coding[key]:
			references.append(key)
			starts.append(i[0])
			ends.append(i[1])

	# #Correct the reference if necessary:
	# for idx, item in enumerate(references):
	# 	if item == "Pf_M76611":
	# 		references[idx] = "M76611"
	# 	elif item == "Pf3D7_API_v3":
	# 		references[idx] = "PFC10_API_IRAB"
	# 	else:
	# 		pass
		

	coverage = []
	for index in tqdm(range(len(references))):
		cov = sam.count_coverage(reference=references[index], start=int(starts[index]), end=int(ends[index]))
		total_cov = numpy.sum(cov, axis=0)
		coverage.append(total_cov) ## Assuming coding regions do NOT OVERLAP

	coverage = numpy.concatenate(coverage, axis=0)
	with open("Coverage_Results.txt", "a") as results_file:
		results_file.write("File :{}\n\n" .format(bamfile) +
							"Coding coverage: {}\n" .format(coverage) +
							"Coding mean coverage: {}\n\n" .format(numpy.mean(coverage))) 
	print "Coding coverage: {}" .format(coverage)
	print "Coding mean coverage: {}" .format(numpy.mean(coverage))

	nc_references = []
	nc_starts = []
	nc_ends = []
	nc_coverage = []

	for key in noncoding:
		for i in noncoding[key]:
			nc_references.append(key)
			nc_starts.append(i[0])
			nc_ends.append(i[1])


	for index in tqdm(range(len(nc_references))):
		nc_cov = sam.count_coverage(reference=nc_references[index], start=int(nc_starts[index]), end=int(nc_ends[index]))
		nc_total_cov = numpy.sum(nc_cov, axis=0)
		nc_coverage.append(nc_total_cov)

	nc_coverage = numpy.concatenate(nc_coverage, axis=0)
	with open("Coverage_Results.txt", "a") as results_file:
		results_file.write("Non-Coding coverage: {}\n" .format(nc_coverage) +
							"Non-Coding mean coverage: {}\n\n" .format(numpy.mean(nc_coverage)))
	print "Non-Coding coverage: {}" .format(nc_coverage)
	print "Non-Coding mean coverage: {}" .format(numpy.mean(nc_coverage))

# Make funtion self-executable:

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		#get_coding_coverage(element)
		#get_flags(element)
		#get_unaligned_mate(element)
		#get_unaligned(element)
		get_frag_len(element)
		get_MAPQ(element)