#!/usr/bin/env python

import csv
import sys

def from_anno_to_table(csv_file):

	with open(str(csv_file).replace("hits.csv", "table.csv"), "a+") as file:
		file.write("Gene\tGene-Coverage\tGene-size\tPre_Cov\tPre_size\tPost_Cov\tPost_size\tChrom\n")

	with open(csv_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if not row: continue
			if row[0].startswith("Pf3D7"):
				chrom = row[0]
				pass
			elif row[0].startswith("ID="):
				 with open(str(csv_file).replace("hits.csv", "table.csv"), "a+") as file:
				 	for element in row:
				 		if "ID=" in element:
				 			file.write(element.replace("ID=","")+"\t")
				 		else:
				 			file.write(element+"\t")
				 	file.write(chrom+"\n")


if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		from_anno_to_table(element)
