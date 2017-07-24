#!/usr/bin/env python

import sys
import pysam
import numpy
import csv
import pandas as pd

depths = {"10G":4535787, "1.2B":4997472, "A7K9": 10915415, "C2": 5176848, "E5K9": 9366386}

def get_coverage(bamfile, chrom, start, stop):
	bam = pysam.AlignmentFile(bamfile, "rb")
	cov = bam.count_coverage(reference=chrom, start=start, end=stop)
	total_cov_vect = numpy.sum(cov, axis=0)
	total_cov = sum(total_cov_vect)

	return total_cov
	

def check_peaks(csv_file):
	
	csv_file = csv_file.replace("./", "")
	bam1 = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/{}_me_sort_q5.bam" .format(csv_file.split("_")[0])
	bam2 = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/{}_me_sort_q5.bam" .format(csv_file.split("_")[1])
	
	depth_factor = float(depths[csv_file.split("_")[0]])/depths[csv_file.split("_")[1]]

	with open(csv_file.replace(".csv", ".bed"), 'rb') as csvfile:
		peaks = csv.reader(csvfile, delimiter='\t')
		calc_fe = []
		for line in peaks:		
			cov1 = get_coverage(bam1, line[0], int(line[1]), int(line[2]))
			cov2 = get_coverage(bam2, line[0], int(line[1]), int(line[2]))
			calc_fe.append(float(cov1)/(cov2*depth_factor))

	input_file = pd.read_csv(csv_file, sep='\t')
	input_file["Calculated_FE"] = calc_fe
	input_file.to_csv(csv_file.replace(".csv", "_calcFE.csv"), index=False, sep='\t')

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		check_peaks(element)

