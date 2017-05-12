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
chrom2 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_02_v3"]
chrom3 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_03_v3"]
chrom4 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_04_v3"]
chrom5 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_05_v3"]
chrom6 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_06_v3"]
chrom7 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_07_v3"]
chrom8 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_08_v3"]
chrom9 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_09_v3"]
chrom10 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_10_v3"]
chrom11 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_11_v3"]
chrom12 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_12_v3"]
chrom13 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_13_v3"]
chrom14 = sorted_gff.loc[sorted_gff[0] == "Pf3D7_14_v3"]
chromAPI = sorted_gff.loc[sorted_gff[0] == "Pf3D7_API_v3"]
chromM76611 = sorted_gff.loc[sorted_gff[0] == "Pf_M76611"]


chromlist = [chrom1, chrom2, chrom3, chrom4, chrom5, chrom6, chrom7, chrom8, chrom9, chrom10, chrom11, chrom12, chrom13, chrom14, chromAPI, chromM76611]

for chrom in chromlist:
	
	chrom1_genes = chrom.loc[chrom[2] == "gene"]

	chrom1_refvect = pandas.DataFrame(range(1, int(sizes[chrom[0].iloc[0]])))
	empty = [0]*(int(sizes[chrom[0].iloc[0]])-1)
	chrom1_refvect[1] = empty
	chrom1_refvect[2] = empty

	for index, row in chrom1_genes.iterrows():
		chrom1_refvect.loc[row[3]:row[4],1] = row[8].split(";")[1][12:]

	for index,row in bed.iterrows():
		chrom1_refvect.loc[row[1]:row[2],2] = 1

	gen_cov = chrom1_refvect.loc[:,1:2].groupby([1]).sum()

	gen_total = chrom1_refvect[1].value_counts()

	result = gen_cov[2]/gen_total*100

	print "\n{}\n" .format(chrom[0].iloc[0])
	print result[result != 0]
	print "\n------------------------------------------------------------------------------------------"



