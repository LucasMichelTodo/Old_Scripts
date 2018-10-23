#!/usr/bin/env python

import sys
import os
import subprocess

def align(filenames):

	subprocess.call("mkdir Results", shell=True)

	tupple_list = []
	i = 0
	while i < len(filenames):
		tupple_list.append((filenames[i], filenames[i+1]))
		i += 2

	for tupple in tupple_list:

		# cmd = "bbduk.sh in={} in2={} out1={} out2={} outm={} ref=/home/lucas/Programs/bbmap/resources/adapters.fa  ktrim=r k=22 mink=6 overwrite=t"\
		#  	.format(tupple[0], tupple[1], tupple[0].replace(".fastq", "_clean.fastq"), tupple[1].replace(".fastq", "_clean.fastq"), tupple[0].split("_")[0]+"_bad.fastq")
		# subprocess.call(cmd, shell=True)

		cmd = "~/Programs/bowtie2-2.3.0-legacy/bowtie2 -p 4 -t --very-sensitive --local -5 4 -3 4 -I 50 -X 200 -x ~/Programs/bowtie2-2.3.0-legacy/Pf3D7 -1 {} -2 {} > ./Results/{}"\
			.format(tupple[0],
					tupple[1],
					tupple[0].replace("_read1_clean.fastq.gz", ".sam"))

		subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    filenames = sys.argv[1:]
    align(filenames)
