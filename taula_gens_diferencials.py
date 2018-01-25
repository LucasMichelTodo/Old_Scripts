#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rosetta_to_dict import *

with open("/home/lucas/ISGlobal/Gen_Referencies/Gens_diferencials_taula.txt", "w+") as outfile:
	outfile.write("ID\tOld_Name\tOld_Family\tExp\tExp_in10G\tVariant\tClassificació\tAnnotation\n")

with open("/home/lucas/ISGlobal/Gen_Referencies/Taula_gens_diferencials.txt", "r+") as file2:
	is_1st_line = True
	for line in file2:
		if is_1st_line:
			is_1st_line = False
		else:
			gene = line.split()[0]
			if line.split()[3] == "FALSE":
				clasif = "Regular"
			elif line.split()[3] == "TRUE" and line.split()[2] == "TRUE":
				clasif = "Variant-active"
			elif line.split()[3] == "TRUE" and line.split()[2] == "FALSE":
				clasif = "Variant-inactive"
			print gene
			with open("/home/lucas/ISGlobal/Gen_Referencies/Gens_diferencials_taula.txt", "a+") as outfile:
				outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" .format(gene, rosetta[gene]["old_names"], 
													rosetta[gene]["old_fam"], 
													line.split()[1], 
													line.split()[2], 
													line.split()[3],
													clasif,
													rosetta[gene]["annot"]))

