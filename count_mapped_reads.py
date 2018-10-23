#!/usr/bin/env python

import sys
import subprocess

nreads = {}

def count_mapped_reads(bamfile):

	# # Command for single-end reads
	# cmd = "samtools view -F 0x904 -c {}"


	# Command for paired-end reads
	 #""| cut -f 1 | sort | uniq | wc -l"


	cmd = "samtools view -F 0x4 {} " .format(bamfile)
	cmd = cmd.split()

	output1 = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]

	cmd =  "cut -f 1"
	cmd = cmd.split()
	output2 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(output1)[0]

	cmd =  ["sort"]
	output3 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(output2)[0]

	cmd =  ["uniq"]
	output4 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(output3)[0]

	cmd =  "wc -l"
	cmd = cmd.split()
	output5 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(output4)[0]

	nreads[bamfile] = output5.strip()


if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		count_mapped_reads(file)

for bam, counts in nreads.iteritems():
	print bam, counts
