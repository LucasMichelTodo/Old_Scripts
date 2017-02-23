#!/usr/bin/env python

import sys
import os
import subprocess
from tqdm import tqdm

filenames = sys.argv[1:]

for file in tqdm(filenames):
	print file
	if file[-4:] == ".sam":
		cmd = "samtools stats {} > {}" .format(file, file.replace(".sam", "_stats.txt"))
		subprocess.call(cmd, shell=True)
		print cmd
	elif file[-4:] == ".bam": 
		cmd = "samtools stats {} > {}" .format(file, file.replace(".bam", "_stats.txt"))
		subprocess.call(cmd, shell=True)
	else:
		print "Not a bam or sam file!"