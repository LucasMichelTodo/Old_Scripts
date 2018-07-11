#!/usr/bin/env python

# Import packages 
import pandas
import csv
import sys
from tqdm import tqdm

# Create a dictionary with chromosome lengths
sizes = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
		else:
			pass

#Import GFF file for annotation.
gff = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", sep="\t", header=None, index_col=None, skiprows=18)

#Sort GFF 
sorted_gff = gff.sort_values([0,3])

#Separate chromosomes.
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

#Create a chromosome list to iterate over.
chromlist = [chrom1, chrom2, chrom3, chrom4, chrom5, chrom6, chrom7, chrom8, chrom9, chrom10, chrom11, chrom12, chrom13, chrom14, chromAPI, chromM76611]

#Main function
def annotate_bed(bed_file):

	#Load bedfile
	bed = pandas.read_csv(filepath_or_buffer= bed_file, sep="\t", header=None, index_col=None, skiprows=1)
	#bed = pandas.read_csv(filepath_or_buffer=csv_file, sep="\t", header=0, index_col=0)

	for chrom in chromlist:

		#Select only fields in the gff that correspond to genes.	
		chrom_genes = chrom.loc[chrom[2] == "gene"]

		#Create a df of the length of the chromosome and append two empty columns.
		chrom_vect = pandas.DataFrame(range(1, int(sizes[chrom[0].iloc[0]]))) 
		empty = [0]*(int(sizes[chrom[0].iloc[0]])-1)
		chrom_vect[1] = empty
		chrom_vect[2] = empty

		#Span gene names over every position they occupy in the genome. (on empty vector 1)
		for index, row in chrom_genes.iterrows(): 
			chrom_vect.iloc[row[3]:row[4],1] = row[8].split(";")[0:1]

		#Span 1s over every covered region in the genome. (on empty vector 2)
		for index,row in bed.iterrows():
			if row[0] == chrom[0].iloc[0]:
				chrom_vect.iloc[row[1]:row[2],2] = 1

		#Count number of 1s per gene.
		gen_cov = chrom_vect.loc[:,1:2].groupby([1]).sum()

		#Count gene length.
		gen_total = chrom_vect[1].value_counts()

		#Divide 1s per gene / gene lengths to get coverage percentage.
		result = gen_cov[2]/gen_total*100

		#Look look at pre and post regions of every gene: (1000bp/previous or next gene)
		surround = []

		for element in tqdm(result.index):

			start = chrom_vect.loc[chrom_vect[1] == element].index[0]
			stop = chrom_vect.loc[chrom_vect[1] == element].index[-1]
			i = 1
			pre_cov = 0
			preregion = 0
			postregion = 0
			post_cov = 0

			# Count coverage in previous 1000bp/gene:
			while i < 1001:
				if chrom_vect.iloc[start-i,1] == 0:
					if chrom_vect.iloc[start-i,2] == 1:
						pre_cov += 1
						preregion += 1
						i += 1
					else:
						preregion += 1
						i += 1
				else:
					i = 1001


			# Count converage in following 1000bp/gene:
			i = 1

			while i < 1001:
				if stop+i < int(sizes[chrom[0].iloc[0]])-1:
					
					if chrom_vect.iloc[stop+i,1] == 0:
						if chrom_vect.iloc[stop+i,2] == 1:
							post_cov += 1
							postregion += 1
							i += 1
						else:
							postregion += 1
							i += 1
					else:
						i = 1001
				else:
					i = 1001

			surround.append([element,[pre_cov, preregion],[post_cov,postregion]])
 

		#Create a df summarizing all results obtained:

		a = [x[0] for x in [y[1] for y in surround]]
		b = [x[1] for x in [y[1] for y in surround]]
		c = [x[0] for x in [y[2] for y in surround]]
		d = [x[1] for x in [y[2] for y in surround]]

		results_df = pandas.DataFrame({'Gene-Coverage':gen_cov[2].tolist(), 'Gene-size':gen_total, 'Pre_Cov':a, 'Pre_size':b, 'Post_Cov':c, 'Post_size':d}, columns=['Gene-Coverage','Gene-size','Pre_Cov','Pre_size','Post_Cov','Post_size'], index=result.index)

		#Print results_df into a csv:
		with open(bed_file.replace(bed_file[-4:],"_annotation.csv"), "a+") as csv_file:
			csv_file.write("\n{}\n" .format(chrom[0].iloc[0]))
		results_df.to_csv(path_or_buf=bed_file.replace(bed_file[-4:],"_annotation.csv"), sep="\t", mode="a+")

		#Print only rows with a non-zero coverage somwhere and eliminate first row:
		no_zero_results = results_df[(results_df['Gene-Coverage'] != 0) | (results_df['Pre_Cov'] != 0) | (results_df['Post_Cov'] != 0)]
		no_zero_results = no_zero_results.iloc[1:]

		#Print positive results to a csv:
		with open(bed_file.replace(bed_file[-4:],"_annotation_hits.csv"), "a+") as csv_file:
			csv_file.write("\n{}\n" .format(chrom[0].iloc[0]))
		no_zero_results.to_csv(path_or_buf=bed_file.replace(bed_file[-4:],"_annotation_hits.csv"), sep="\t", mode="a+")


#If called from command line, execute main function for each file.
if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		annotate_bed(element)