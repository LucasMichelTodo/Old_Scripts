#!/usr/bin/env python

import sys
import subprocess

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames

for file in filenames:
	cmd = "samtools view -b -q 5 {} > {}" .format(file, file.replace(".bam", "_q5.bam"))
	print cmd 
	subprocess.call(cmd, shell=True)



