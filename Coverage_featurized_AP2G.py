#!/usr/bin/env python

import sys
import pysam
import numpy
numpy.set_printoptions(threshold='nan')
import csv
import pandas as pd

#depths = {"10G":10375962.0, "1.2B":12066703.0, "A7K9": 30393549.0, "C2": 12006887.0, "E5K9": 25195030.0} #A usar per metilacio
#depths = {"10G":30798124.0, "1.2B":30009498.0, "A7K9": 17489454.0, "C2": 32059287.0, "E5K9": 21544493.0} #A usar per inputs
depths = {"10G":22659880.0, "1.2B":30648669.0, "C2": 21165995.0} #A usar per acetilacio

def get_coverage(bamfile, chrom, start, stop):
	bam = pysam.AlignmentFile(bamfile, "rb")
	cov = bam.count_coverage(reference=chrom, start=start, end=stop)
	total_cov_vect = numpy.sum(cov, axis=0)
	total_cov = sum(total_cov_vect) / float(stop-start)

	return total_cov, total_cov_vect
	

def featurized_coverage(bamfile):
	
	with open(bamfile.replace("_sort_q5.bam", "_cov_AP2G.csv"), "w+") as file2: #Ensure we create a new file "from scratch" every time the program is run.
		file2.write("Gene\tGene-cov\t5-cov\t3-cov\tChrom\tAnnotations\n")

	with open("/home/lucas/ISGlobal/Gen_Referencies/AP2G.gff", "r") as file:  
		gff_file = csv.reader(file, delimiter="\t")
		
		for line in gff_file:
			pre = get_coverage(bamfile, line[0], int(line[3]), int(line[3])+int(line[9]))[0] / depths[bamfile.split("_")[0]]*1000000.0
			post = get_coverage(bamfile, line[0], int(line[4])-int(line[10]), int(line[4]))[0] / depths[bamfile.split("_")[0]]*1000000.0
			covGene = get_coverage(bamfile, line[0], int(line[3])+int(line[9]), int(line[4])-int(line[10]))[0] / depths[bamfile.split("_")[0]]*1000000.0

			if line[6] == "+":
				cov5 = pre
				cov3 = post
			else:
				cov5 = post
				cov3 = pre

			with open(bamfile.replace("_sort_q5.bam", "_cov_AP2G.csv"), "a+") as file2:
				file2.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(covGene)+"\t"+str(cov5)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n")

			
			with open(bamfile.replace("_sort_q5.bam", "_cov_vector_pre_AP2G.csv"), "a+") as file3:
				vect = (get_coverage(bamfile, line[0], int(line[3]), int(line[3])+int(line[9]))[1] / depths[bamfile.split("_")[0]])*1000000.0
				for element in vect:
					file3.write(str(element)+",")

			with open(bamfile.replace("_sort_q5.bam", "_cov_vector_post_AP2G.csv"), "a+") as file3:
				vect = (get_coverage(bamfile, line[0], int(line[4])-int(line[10]), int(line[4]))[1] / depths[bamfile.split("_")[0]])*1000000.0
				for element in vect:
					file3.write(str(element)+",")

			with open(bamfile.replace("_sort_q5.bam", "_cov_vector_ORF_AP2G.csv"), "a+") as file3:
				vect = (get_coverage(bamfile, line[0], int(line[3])+int(line[9]), int(line[4])-int(line[10]))[1] / depths[bamfile.split("_")[0]])*1000000.0
				for element in vect:
					file3.write(str(element)+",")



if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		featurized_coverage(element)