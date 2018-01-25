#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import subprocess as sp

### Load gff and elongate each gene 1000bp pre and post.
# # Load GFF for annotation
# ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff") 

# # Sort GFF
# ref = ref.sort()

# # Filter only those annotation that correspond to genes.
# gene_ref = ref.filter(lambda x: x[2] == "gene")

# # Add 1kb before and after each gene to capture 5' and 3' regions.
# gene_ref = gene_ref.slop(b=1000, g="/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome")


### Load directly edited gff (by elongator.r) with genes elongates 1000bp or until next/previous gene.
gene_ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes.gff")

# Function for calculating overlap.
def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

# Anotating main function.
def annotate_from_genes(bed_file):

	with open(bed_file.replace(".bed", "_annotated.csv"), "w+") as file: # Creating header for annotation file. Setting "w+" ensures the file will be created from scratch if re-run.
		file.write("Gene\tGene-cov\t5-cov\t3-cov\tChrom\tAnnotations\n")

	bed = py.BedTool(bed_file) # Load bed file with peaks. 
	intersect = gene_ref.intersect(bed, wo=True) # Intersect each gene in the gff with the peaks in the bed file.

	for line in intersect: # Calculate overlap between first 1kb, gene ORF and last 1kb with the intersecting peak.
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
			file.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(gene)+"\t"+str(cov5)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n") # Write results into annotation file.


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		annotate_from_genes(element)