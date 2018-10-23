#!/usr/bin/env python

import sys
import pysam
import numpy
numpy.set_printoptions(threshold='nan')
import csv
import pandas as pd



# Load file with number of reads into a dictionary:
nreads = {}
#with open("/home/lucas/ISGlobal/Chip_Seq/Noves_dades/chip_nou_nreads.txt", "r+") as infile:  ## New Chip
with open("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/chip_vell_nreads.txt", "r+") as infile: ## Cristina's Old Chip
	for line in infile:
		nreads[line.split()[0]] = float(line.strip().split()[1])

def get_coverage(bamfile, chrom, start, stop):
	bam = pysam.AlignmentFile(bamfile, "rb")
	cov = bam.count_coverage(reference=chrom, start=start, end=stop)
	total_cov_vect = numpy.sum(cov, axis=0)
	total_cov = sum(total_cov_vect) / float(stop-start)

	return total_cov, total_cov_vect

# Function for dividing a list into "n" roughly equal length lists:
def split(a, n):
	k, m = divmod(len(a), n)
	splits = [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)]
	return [[i[0], i[-1]] for i in splits]


def featurized_coverage(bamfile):

	with open(bamfile.replace("_sort_q5.bam", "_cov_intervals.csv"), "w+") as file2: #Ensure we create a new file "from scratch" every time the program is run.
		file2.write("Gene\tSeg1\tSeg2\tSeg3\tSeg4\tSeg5\tSeg6\tSeg7\tSeg8\tChrom\tAnnotations\n")

	with open("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff", "r") as file:
		gff_file = csv.reader(file, delimiter="\t")

		for line in gff_file:
			gene_range = range(int(line[3]), int(line[4]))
			segments = split(gene_range, 8)
			seg_cov = []

			for seg in segments:
				seg_cov.append(get_coverage(bamfile, line[0], seg[0], seg[1])[0] / nreads[bamfile]*1000000.0)

			if line[6] == "-":
				seg_cov.reverse()

			with open(bamfile.replace("_sort_q5.bam", "_cov_intervals.csv"), "a+") as file2:
				coverages = "\t".join(str(i) for i in seg_cov)
				file2.write(line[8].split(";")[0].replace("ID=","")+"\t"+coverages+"\t"+line[0]+"\t"+line[8]+"\n")

			print line[8].split(";")[0].replace("ID=","")
			print "----------------"



if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		featurized_coverage(element)
