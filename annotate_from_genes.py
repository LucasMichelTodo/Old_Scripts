#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import subprocess as sp

ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff")
ref = ref.sort()
gene_ref = ref.filter(lambda x: x[2] == "gene")
gene_ref = gene_ref.slop(b=1000, g="/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome")

def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

def annotate_from_genes(bed_file):

	with open(bed_file.replace(".bed", "_annotated.csv"), "w+") as file:
		file.write("Gene\tGene-Coverage\t5-cov\t3-cov\tChrom\tAnnotations\n")

	bed = py.BedTool(bed_file)
	intersect = gene_ref.intersect(bed, wo=True)

	for line in intersect:
		pre = round(overlap(int(line[3]), int(line[3])+1000, int(line[10]), int(line[11])) / float(10), 2)
		gene = round(overlap(int(line[3])+1000, int(line[4])-1000, int(line[10]), int(line[11])) / float((int(line[4])-1000) - (int(line[3])+1000))*100, 2)
		post = round(overlap(int(line[4])-1000, int(line[4]), int(line[10]), int(line[11])) / float(10), 2)
		
		if line[6] == "+":
			cov5 = pre
			cov3 = post
		else:
			cov5 = post
			cov3 = pre

		with open(bed_file.replace(".bed", "_annotated.csv"), "a+") as file:
			file.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(gene)+"\t"+str(cov5)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n")


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		annotate_from_genes(element)