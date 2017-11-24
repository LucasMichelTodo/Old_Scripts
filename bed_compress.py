#!/usr/bin/env python

import sys

def compress_bed(input_bed):
	with open(input_bed, "r+") as infile:
		i = 0
		pos = 0
		cov = []
		for line in infile:
			if i == 0:
				chrom_ref = line.split()[0]
				cov.append(line.split()[3])
				i += 1
				pos += 1
				
			if i !=0 and i < 5:
				chrom = line.split()[0]
				if chrom == chrom_ref:
					cov.append(line.split()[3])
					i += 1
					pos += 1
				else:
					i = 1
					with open(input_bed.replace(".txt", "_compressed.bed"), "a+") as oufile:
						outfile.write(chrom+"\t"+"")
					cov = line.split()[3]


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		compress_bed(element)