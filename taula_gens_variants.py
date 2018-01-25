#!/usr/bin/env python

from rosetta_to_dict import *

gens_variants = []

with open("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants.txt", "r+") as file2:
	for line in file2:
		for key, value in rosetta.iteritems():
			if line.split("\t")[0].strip() in value["old_names"]:
				gens_variants.append(key)

with open("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff", "r+") as file3:
	for line in file3:
		if line.split("\t")[8].split(";")[0].replace("ID=", "") in gens_variants:
			print line.split("\t")[8].split(";")[0].replace("ID=", "") + "---->" +line.split("\t")[8].split(";")[1].replace("description=", "")
			with open("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_gff2.txt", "a+") as outfile:
				outfile.write(line.split("\t")[8].split(";")[0].replace("ID=", "")+"\t"+line)