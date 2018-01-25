#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import subprocess as sp
from Coverage_featurized2 import get_coverage

gene_ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff")

# Anotating main function.
def annotate_from_genes(bamfile):

	with open(bamfile.replace(".bam", "_annotated_cov.csv"), "w+") as file: # Creating header for annotation file. Setting "w+" ensures the file will be created from scratch if re-run.
		file.write("Gene\t5_cov\tORF_cov\t3_cov\tChrom\tAnnotations\n")

	for line in gene_ref: 
		
		print line
		
		pre = get_coverage(bamfile, line[0], int(line[3]), int(line[3])+int(line[9]))
		ORF = get_coverage(bamfile, line[0], int(line[3])+int(line[9]), int(line[4])-int(line[10]))
		post = get_coverage(bamfile, line[0], int(line[4])-int(line[10]), int(line[4]))

		if line[6] == "+":
			cov5 = pre
			cov3 = post
		else:
			cov5 = post
			cov3 = pre
		

		with open(bamfile.replace(".bam", "_annotated_cov.csv"), "a+") as file:
			file.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(cov5)+"\t"+str(ORF)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n")


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for file in filenames:
		annotate_from_genes(file)