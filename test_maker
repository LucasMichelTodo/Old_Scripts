#!/usr/bin/env python

import sys
import os
import subprocess
from tqdm import tqdm

def samplify(nsample, rawfile):
	if rawfile[-4:] == ".sam":
		cmd = "head -{} {} > {}" .format(nsample, rawfile, rawfile.replace(".sam", "_sample.sam"))
		subprocess.call(cmd, shell = True)
	elif rawfile[-4:] == ".bam": 
		cmd = "head -{} {} > {}" .format(nsample, rawfile, rawfile.replace(".bam", "_sample.bam"))
		subprocess.call(cmd, shell = True)
	elif rawfile[-6:] == ".fastq": 
		cmd = "head -{} {} > {}" .format(nsample, rawfile, rawfile.replace(".fastq", "_sample.fastq"))
		subprocess.call(cmd, shell = True)
	else:
		print "Not a bam/sam or fastq file!"

sample_size = sys.argv[1]
filenames = sys.argv[2:]

for file in tqdm(filenames):
	samplify(sample_size, file)

# USAGE
#
# Call "test_maker <int> <file/s>" it will create a file with only first <int> lines for each <file>