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

# Get reference annotated genome

ref_fasta = open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-30_Pfalciparum3D7_AnnotatedCDSs.fasta", "r+")
coding_regions = {}
for line in ref_fasta:
	if line.startswith(">"):
		gene = line.split("|")
		gene = map(lambda x: x.strip(" >\n"), gene)
		dict_string = map(lambda x: x.split("="), gene[1:])
		gene_dict = {}
		for element in dict_string:
			gene_dict[element[0]] = element[1]
		coding_regions[gene[0]] = gene_dict
coding_df = pandas.DataFrame.from_dict(data=coding_regions, orient = "index")
#print df

def get_frag_len(samfile):
	sam = pysam.AlignmentFile(samfile, "r")
	lens = []
	for read in sam.fetch():
		lens.append(read.template_length)
	with open(samfile.replace(".sam", "_lengths.csv"),"w") as csvfile:
		writer = csv.writer(csvfile, delimiter = "\t", quoting=csv.QUOTE_MINIMAL)
		writer.writerow(lens)

def get_MAPQ(samfile):
	sam = pysam.AlignmentFile(samfile, "r")
	mapq = []
	for read in sam.fetch():
		mapq.append(read.mapping_quality)
	with open(samfile.replace(".sam", "_MAPQ.csv"), "w") as csvfile:
		writer = csv.writer(csvfile, delimiter = "\t", quoting=csv.QUOTE_MINIMAL)
		writer.writerow(mapq)

def get_unaligned(samfile):
	cmd = "samtools view -h -f 4 {} > {}" .format(samfile, samfile.replace(".sam", "_unmapped.sam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools view -h -F4 -f8 {} > {}" .format(samfile, samfile.replace(".sam", "_mateunmapped.sam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools merge -h {} {} {} {}" .format(samfile, samfile.replace(".sam","_somethingununmapped.sam"), 
												samfile.replace(".sam", "_unmapped.sam"),
												samfile.replace(".sam", "_mateunmapped.sam"))
	subprocess.call(cmd, shell=True)

	
def get_unaligned_mate(samfile):
	sam = pysam.AlignmentFile(samfile, "r")
	unpaired_mates = 0
	for read in sam.fetch():
		if read.mate_is_unmapped:
			unpaired_mates += 1
	print unpaired_mates

def get_flags(samfile):
	flags = {}
	sam = pysam.AlignmentFile(samfile, "r")
	for read in sam.fetch():
		if read.flag in flags:
			flags[read.flag] += 1
		else:
			flags[read.flag] = 1
	print flags

def get_coding_coverage(bamfile):
	sam = pysam.AlignmentFile(bamfile, "rb")
	
	regions = coding_df["location"]
	references = []
	starts = []
	ends = []
	for element in regions:
		references.append(element.split(":")[0])
		starts.append(element.split(":")[1].strip("([+-])").split("-")[0])
		ends.append(element.split(":")[1].strip("([+-])").split("-")[1])

	coverage = [] 
	index = 0
	while index <= len(references)-1:
		cov = sam.count_coverage(reference=references[index], start=int(starts[index]), end=int(ends[index]))
		total_cov = numpy.sum(cov, axis=0)
		coverage.append(total_cov)
		index += 1

	coverage = numpy.concatenate(coverage, axis=0)
	print coverage
	print numpy.mean(coverage)



filenames = sys.argv[1:]
print filenames
for element in tqdm(filenames):
	get_coding_coverage(element)
	#get_flags(element)
	#get_unaligned_mate(element)
	#get_unaligned(element)
	#get_frag_len(element)
	#get_MAPQ(element)