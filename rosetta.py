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

genes = {}

def translate(file):
	with open(file) as file:
		for line in file:
			genes[(line.strip())] = []

		for i in tqdm(genes):
			if i in str(rosetta.values()):
				for key, value in rosetta.iteritems():
					if i in value:
							# print i
							# print value
							# print "------------------------------------------------------------------"
							genes[i].append(key)
			else:
				genes[i].append("NA")

	with open("rosseta_IDs3.txt", "a+") as outfile:
		for key, value in genes.iteritems():
			outfile.write(key+"\t")
			first = True
			for i in value:
				if first: 
					outfile.write(i)
					first = False
				else:
					outfile.write(","+i)

			outfile.write("\n")



if __name__ == "__main__":
	file = sys.argv[1]
	translate(file)