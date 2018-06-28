#!/usr/bin/env python

import sys
from tqdm import tqdm
from rosetta_annotated_dict import rosetta

genes = []

def translate(filein):
	with open(filein) as file:
		for line in file:
			genes.append([line.strip(),["NA"],"NA", "NA"])

		for i in tqdm(genes):
			for gene in rosetta:
				if i[0] in rosetta[gene]["old_refs"]:
					if i[1][0] == "NA":
						i[1] = [gene]
					else:
						i[1].append(gene)
					i[2] = rosetta[gene]["annot"]
					i[3] = rosetta[gene]["name"]


	with open(filein.replace(".txt", "_rosetta_annotated.txt"), "a+") as outfile:
		for i in genes:
			outfile.write(i[0]+"\t")
			first = True
			for y in i[1]:
				if first:
					outfile.write(y)
					first = False
				else:
					outfile.write(","+y)
			outfile.write("\t"+i[2]+"\t"+i[3]+"\n")


if __name__ == "__main__":
	filein = sys.argv[1]
	translate(filein)
