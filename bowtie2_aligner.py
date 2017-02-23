#!/usr/bin/env python

import sys
import os
import subprocess
from tqdm import tqdm

	
# Set in/out paths
inpath = "/home/lucas/ISGlobal/TestSet/align_tests/inputs/"
outpath = "/home/lucas/ISGlobal/TestSet/align_tests/"

# Set params
params = {"params_1":"-p 4 --very-sensitive --local -5 4 -3 4 -I 50 -X 250", 
		"params_2":"-p 4 --very-sensitive --local -5 4 -3 4 -I 50 -X 2000"}

fd_rawlist = os.listdir(inpath)

for key in params:
	print key

names = []
for file in fd_rawlist:
	names.append(file[:-19])

for x in tqdm(names):
	cmd = "~/Programs/bowtie2-2.3.0-legacy/bowtie2 {} -x ~/Programs/bowtie2-2.3.0-legacy/Pf3D7 -1 {} -2 {} > {}" .format(params[2], inpath+x+"_read1_sample.fastq", inpath+x+"_read2_sample.fastq", outpath+x+".sam")
	subprocess.call(cmd, shell=True)





