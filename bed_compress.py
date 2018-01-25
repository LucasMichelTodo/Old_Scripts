#!/usr/bin/env python

import sys
import numpy
from tqdm import tqdm


# Load a bed file and convert every 5 entries into a single one with mean value of them.
def compress_bed(input_bed):
	with open(input_bed, "r+") as infile:
		i = 0 #Loop counter.
		pos = 0 #Counter for genomic position.
		cov = [] #Container for bed values to average.
		
		for line in tqdm(infile):
			if i == 0: #handling first line
				chrom_ref = line.split()[0]
				cov.append(line.split()[3])
				i += 1
				pos += 1
				
			if i !=0 and i < 5: #Main loop
				chrom = line.split()[0]

				if chrom == chrom_ref: #Check that we are not in a "chromosome change"
					cov.append(line.split()[3])
					i += 1
					pos += 1
					
				else: # If we are in a "chromosome change", write line with average buffered values and restart loop.
					print chrom
					with open(input_bed.replace(".txt", "_compressed.bed"), "a+") as outfile:
						outfile.write(chrom_ref+"\t"+str(pos-i)+"\t"+str(pos+1)+"\t"+str(numpy.array(cov).astype(numpy.float).mean())+"\n")
					i = 1
					pos = 0
					cov = [line.split()[3]]
					chrom_ref = chrom
					
			else:
				with open(input_bed.replace(".txt", "_compressed.bed"), "a+") as outfile:
					outfile.write(chrom_ref+"\t"+str(pos-i)+"\t"+str(pos+1)+"\t"+str(numpy.array(cov).astype(numpy.float).mean())+"\n")
				i = 1
				pos += 1
				cov = [line.split()[3]]
				

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		compress_bed(element)