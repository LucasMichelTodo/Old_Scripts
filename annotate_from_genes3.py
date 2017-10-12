#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import subprocess as sp
from Coverage_featurized2 import get_coverage

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
gene_ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff")

# Function for calculating start and stop of overlapping sequence. Returns (0,0) if there is no overlap.
def overlap(min1, max1, min2, max2):
	vect = (max(min1, min2), min(max1, max2))
	if vect[0] > vect[1]:
		vect = (0,0)
	return vect

# Anotating main function.
def annotate_from_genes(bed_file, bamfile):

	with open(bed_file.replace(".bed", "_annotated_cov.csv"), "w+") as file: # Creating header for annotation file. Setting "w+" ensures the file will be created from scratch if re-run.
		file.write("Gene\tGene-cov\t5-cov\t3-cov\tChrom\tAnnotations\n")

	bed = py.BedTool(bed_file) # Load bed file with peaks. In order to cover ALL peaks must be a merged file with peaks coming both from 1.2B and 10G!
	intersect = gene_ref.intersect(bed, wo=True) # Intersect each gene in the gff with the peaks in the bed file.

	for line in intersect: # Calculate overlap between pre-region (1kb or until previous gene), gene ORF and post-region with the intersecting peak.
			
		pre = overlap(int(line[3]), int(line[3])+int(line[9]), int(line[12]), int(line[13]))
		if pre[0]+pre[1] > 0:
			pre_cov = get_coverage(bamfile, line[0], pre[0], pre[1])
		else:
			pre_cov = 0
		
		gene = overlap(int(line[3])+int(line[9]), int(line[4])-int(line[10]), int(line[12]), int(line[13]))
		if gene[0]+gene[1] > 0:
			gene_cov = get_coverage(bamfile, line[0], gene[0], gene[1])
		else:
			gene_cov = 0

		post = overlap(int(line[4])-int(line[10]), int(line[4]), int(line[12]), int(line[13]))
		if post[0]+post[1] > 0:
			post_cov = get_coverage(bamfile, line[0], post[0], post[1])
		else:
			post_cov = 0

		if line[6] == "+":
			cov5 = pre_cov
			cov3 = post_cov
		else:
			cov5 = post_cov
			cov3 = pre_cov

		with open(bed_file.replace(".bed", "_annotated_cov.csv"), "a+") as file:
			file.write(line[8].split(";")[0].replace("ID=","")+"\t"+str(gene_cov)+"\t"+str(cov5)+"\t"+str(cov3)+"\t"+line[0]+"\t"+line[8]+"\n") # Write results into annotation file.


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	annotate_from_genes(filenames[0], filenames[1])