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

for file in (filenames):
	if file[-4:] == ".sam":
		stats_file = file.replace(".sam", "_stats.txt")
	else:
		stats_file = file.replace(".bam", "_stats.txt")
	subprocess.call("grep ^SN {} | cut -f 2- > {}" .format(stats_file, stats_file.replace(".txt", "_summary.csv")), shell=True)