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

# sizes["Pf_M76611"] = sizes.pop("M76611")
# sizes["Pf3D7_API_v3"] = sizes.pop("PFC10_API_IRAB")

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
		noncoding[key].append(int(i[0])-1)
		noncoding[key].append(int(i[1])+1)

noncoding["M76611"] = noncoding.pop("Pf_M76611")
noncoding["PFC10_API_IRAB"] = noncoding.pop("Pf3D7_API_v3")

for key in noncoding:
	noncoding[key].insert(0,1)
	noncoding[key].append(sizes[key])

# Set one 0 to 1 (coding region starts at pos 1):
noncoding["PFC10_API_IRAB"][1] = 1

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

	# Correct the reference:
	for idx, item in enumerate(references):
		if item == "Pf_M76611":
			references[idx] = "M76611"
		elif item == "Pf3D7_API_v3":
			references[idx] = "PFC10_API_IRAB"
		else:
			pass
		

	coverage = []
	for index in tqdm(range(len(references))):
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
		nc_cov = sam.count_coverage(reference=nc_references[index], start=int(nc_starts[index]), end=int(nc_ends[index]))
		nc_total_cov = numpy.sum(nc_cov, axis=0)
		nc_coverage.append(nc_total_cov)


	nc_coverage = numpy.concatenate(nc_coverage, axis=0)
	print "Non-Coding coverage: {}" .format(nc_coverage)
	print "Non-Coding mean coverage: {}" .format(numpy.mean(nc_coverage))


filenames = sys.argv[1:]
print filenames
for element in filenames:
	get_coding_coverage(element)
	#get_flags(element)
	#get_unaligned_mate(element)
	#get_unaligned(element)
	#get_frag_len(element)
	#get_MAPQ(element)