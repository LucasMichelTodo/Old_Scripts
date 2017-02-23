#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
from Bio.Seq import Seq
import subprocess
import pysam
from tqdm import tqdm
import csv


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

filenames = sys.argv[1:]
print filenames
for element in tqdm(filenames):
	#get_flags(element)
	#get_unaligned_mate(element)
	#get_unaligned(element)
	get_frag_len(element)
	#get_MAPQ(element)