#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py
import subprocess as sp

# Load GFF for annotation
ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff")

noncoding = ref.filter(lambda x: x[2] == "ncRNA" or x[2] == "snoRNA" or x[2] == "rRNA" or x[2] == "sn" or x[2] == "tRNA" or x[2] == "three_prime_UTR").saveas("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_noncoding.gff")

bannlist = []
for i in noncoding:
	bannlist.append((i.chrom, i.start))

print len(bannlist)

# Sort GFF
ref = ref.sort().saveas("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Sorted.gff")

# Filter only those annotation that correspond to genes. And write result.
gene_ref = ref.filter(lambda x: x[2] == "gene")

def bann_feature(feature):
	for i in bannlist:
		if feature.chrom == i[0] and feature.start == i[1]:
			feature.start = -2
	return feature

gene_ref_banned = gene_ref.each(bann_feature)

gene_ref_final = gene_ref_banned.filter(lambda x: x[3] != "-1").saveas("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Sorted_filtered.gff")
