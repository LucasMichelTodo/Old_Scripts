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


## Creating non-coding regions reference:

regions = coding_df["location"]

sizes = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
		else:
			pass

sizes["Pf_M76611"] = sizes.pop("M76611")
sizes["Pf3D7_API_v3"] = sizes.pop("PFC10_API_IRAB")

chrom = map(lambda x: x.split(":"), regions)
chroms = {}
for pair in chrom:
	if pair[0] in chroms:
		chroms[pair[0]].append(pair[1])
	else:
		chroms[pair[0]] = [pair[1]]

for key in chroms:
	chroms[key] = map(lambda x: x.strip("([+-])").split("-"), chroms[key])

noncoding = {}
for key in chroms:
	noncoding[key] = []
	for i in chroms[key]:
#if i > noncoding[key][-1]:
		noncoding[key].append(int(i[0])-1)
		noncoding[key].append(int(i[1])+1)




for key in noncoding:
	noncoding[key].insert(0,1)
	noncoding[key].append(sizes[key])

## Remove overlapping genes:

n_iter = [0,1]
for i in n_iter:
	for key in noncoding:
		for i in range(len(noncoding[key])):
			if len(noncoding[key][i-2:i+2]) == 4:
				if noncoding[key][i] < noncoding[key][i-1]:
					# print "Overlap in Chrom: {} position {}" .format(key, i)
					# print noncoding[key][i-2:i+2]
					noncoding[key][i-2:i+2] = [min(noncoding[key][i-2:i+2]),max(noncoding[key][i-2:i+2])]
			else:
				pass

## Check no overlapping genes are still present:

for key in noncoding:
	for i in range(len(noncoding[key])):
		if noncoding[key][i] < noncoding[key][i-1] and i != 0 and i !=1:
			print "Overlap in Chrom: {} position {}" .format(key, i)
			print noncoding[key][i-2:i+2]




### FUNCTIONS

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
	for index in tqdm(range(len(references))):
		#print "{} {} {}" .format(references[index], starts[index], ends[index])
		cov = sam.count_coverage(reference=references[index], start=int(starts[index]), end=int(ends[index]))
		total_cov = numpy.sum(cov, axis=0)
		coverage.append(total_cov) ## Assuming coding regions do NOT OVERLAP

	coverage = numpy.concatenate(coverage, axis=0)
	print "Coding coverage: {}" .format(coverage)
	print "Coding mean coverage: {}" .format(numpy.mean(coverage))

	nc_references = []
	nc_starts = []
	nc_ends = []
	nc_coverage = []
	i = 0
	for key in noncoding:
		for x in noncoding[key]:
			if i == 0:
				nc_references.append(key)
				nc_starts.append(x)
				i = 1
			else:
				nc_ends.append(x)
				i = 0


	for index in tqdm(range(len(nc_references))):
		print "{} {} {}" .format(nc_references[index], nc_starts[index], nc_ends[index])
		nc_cov = sam.count_coverage(reference=nc_references[index], start=int(nc_starts[index]), end=int(nc_ends[index]))
		#print (nc_cov)
		nc_total_cov = numpy.sum(nc_cov, axis=0)
		#print nc_total_cov
		nc_coverage.append(nc_total_cov)
		#print (nc_coverage) 

	nc_coverage = numpy.concatenate(nc_coverage, axis=0)
	print "Non-Coding coverage: {}" .format(nc_coverage)
	print "Non-Coding mean coverage: {}" .format(numpy.mean(nc_coverage))


# nc_references = []
# nc_starts = []
# nc_ends = []
# nc_coverage = []
# i = 0
# for key in noncoding:
# 	for x in noncoding[key]:
# 		if i == 0:
# 			nc_references.append(key)
# 			nc_starts.append(x)
# 			i = 1
# 		if i == 1:
# 			nc_ends.append(x)
# 			i = 0

# print len(nc_references)
# print len(nc_starts)
# print len(nc_ends)

filenames = sys.argv[1:]
print filenames
#for element in filenames:
	#get_coding_coverage(element)
	#get_flags(element)
	#get_unaligned_mate(element)
	#get_unaligned(element)
	#get_frag_len(element)
	#get_MAPQ(element)