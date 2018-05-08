#!/usr/bin/env python

import sys
from tqdm import tqdm

rosetta = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", "r+") as file1:
	for line in file1:
		rosetta[str(line.split("\t")[0].strip())] = line.split("\t")[1:]
	for key, value in rosetta.iteritems():
		value[-1] = value[-1].strip()

# print len(rosetta)
# print rosetta

genes = []

def translate(filein):
	with open(filein) as file:
		for line in file:
			genes.append([line.strip(),""])

		for i in tqdm(genes):
			if i[0] in str(rosetta.values()):
				for key, value in rosetta.iteritems():
					if i[0] in value:
						i[1] = key
			else:
				i[1] = "NA"


	with open(filein.replace(".txt", "_rosetta.txt"), "a+") as outfile:
		for i in genes:
			outfile.write(i[0]+"\t")
			first = True
			for y in i[1:]:
				if first:
					outfile.write(y)
					first = False
				else:
					outfile.write(","+y)

			outfile.write("\n")



if __name__ == "__main__":
	filein = sys.argv[1]
	translate(filein)
