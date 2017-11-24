#!/usr/bin/env python

import sys
import pysam
import numpy
import csv
import pandas as pd

depths = {"10G":4535787, "1.2B":4997472, "A7K9": 10915415, "C2": 5176848, "E5K9": 9366386} # comprovar que corespon a q5

def get_coverage(bamfile, chrom, start, stop):
	bam = pysam.AlignmentFile(bamfile, "rb")
	cov = bam.count_coverage(reference=chrom, start=start, end=stop)
	total_cov_vect = numpy.sum(cov, axis=0)
	total_cov = sum(total_cov_vect) / float(stop-start)

	return total_cov
	

def featurized_coverage(bamfile):
	
	with open(bamfile.replace("mapped_sorted.bam", "cov.csv"), "w+") as file2: #Ensure we create a new file "from scratch" every time the program is run.
		file2.write("Gene\tGene-cov\t5-cov\t3-cov\tChrom\tAnnotations\n")

	with open("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff", "r") as file:  
		gff_file = csv.reader(file, delimiter="\t")
		
		for line in gff_file:
			pre = get_coverage(bamfile, line[0], int(line[3]), int(line[3])+int(line[9])) / depths[bamfile.split("_")[0]]*1000000.0
			post = get_coverage(bamfile, line[0], int(line[4])-int(line[10]), int(line[4])) / depths[bamfile.split("_")[0]]*1000000.0
			covGene = get_coverage(bamfile, line[0], int(line[3])+int(line[9]), int(line[4])-int(line[10])) / depths[bamfile.split("_")[0]]*1000000.0

			if line[6] == "+":
				cov5 = pre
				cov3 = post
			else:
				cov5 = post
				cov3 = pre

			with open(bamfile.replace("mapped_sorted.bam", "cov.csv"), "a+") as file2:
				file2.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(covGene)+"\t"+str(cov5)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n")


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		featurized_coverage(element)