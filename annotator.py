#!/usr/bin/env python

# Import packages 
import pandas
import csv

sizes = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
		else:
			pass

gff = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", sep="\t", header=None, index_col=None, skiprows=18)

bed = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Chip_Seq/DATA/Substractions/Macs_diff/diff_A7_E5_narrowpeaks_c3.0_cond1.bed", sep="\t", header=None, index_col=None, skiprows=1)

sorted_gff = gff.sort_values([0,3])

chrom1 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_01_v3"]

chrom1_genes = chrom1.loc[chrom1[2] == "gene"]

chrom1_refvect = pandas.DataFrame(range(1, int(sizes["Pf3D7_01_v3"])))
empty = [0]*(int(sizes["Pf3D7_01_v3"])-1)
chrom1_refvect[1] = empty
chrom1_refvect[2] = empty

for index, row in chrom1_genes.iterrows():
	chrom1_refvect.loc[row[3]:row[4],1] = row[8].split(";")[1]

for index,row in bed.iterrows():
	chrom1_refvect.loc[row[1]:row[2],2] = 1

gen_cov = chrom1_refvect.loc[:,1:2].groupby([1]).sum()

gen_total = chrom1_refvect[1].value_counts()

print gen_cov[2]/gen_total*100






