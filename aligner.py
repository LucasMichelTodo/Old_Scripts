#!/usr/bin/env python

import sys
import os
import subprocess

subprocess.call("mkdir Results", shell=True)

filenames = sys.argv[1:]

tupple_list = []
i = 0
while i < len(filenames):
	tupple_list.append((filenames[i], filenames[i+1]))
	i += 2

for tupple in tupple_list:
	cmd = "~/Programs/bowtie2-2.3.0-legacy/bowtie2 -p 4 --very-sensitive --local -5 4 -3 4 -I 50 -X 200 -x ~/Programs/bowtie2-2.3.0-legacy/Pf3D7 -1 {} -2 {} > ./Results/{}"\
		.format(tupple[0], tupple[1], tupple[0].replace("_read1.fastq", ".sam"))
	subprocess.call(cmd, shell=True)

