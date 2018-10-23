#!/usr/bin/env python

import sys
from tqdm import tqdm
from rosetta_annotated_dict import rosetta

genes = []

def translate(filein):
	with open(filein) as file:
		for line in file:
			genes.append(line.strip().split(":")[0].replace("/", "_"))

		for i in tqdm(genes):
			try:
				anot = rosetta[i]["annot"]
				name = rosetta[i]["name"]
				with open(filein.replace(".txt", "_rosetta_annotated.txt"), "a+") as outfile:
						outfile.write(i+"\t"+anot+"\t"+name+"\n")
			except:
				anot = "NA"
				name = "NA"
				print "Gene {} was not found in database!" .format(i)
				with open(filein.replace(".txt", "_rosetta_annotated.txt"), "a+") as outfile:
						outfile.write(i+"\t"+anot+"\t"+name+"\n")



if __name__ == "__main__":
	filein = sys.argv[1]
	translate(filein)
