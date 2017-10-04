#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import pandas as pd

# Load GFF for annotation
ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff")
ref = ref.slop(l=1000, r=0, g="/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome", s=True)

# Sort GFF
ref = ref.sort()

# Filter only those annotation that correspond to genes.
gene_ref = ref.filter(lambda x: x[2] == "gene")

def get_gene_coverage(bam_file):
	bam = py.BedTool(bam_file)
	bed = bam.bam_to_bed()
	result = gene_ref.coverage(bed, mean=True)
	with open(bam_file.replace(".bam", "_CovMtx_5ext.csv"), "w+") as result_file:
		result_file.write(str(result))

	result

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		get_gene_coverage(element)
