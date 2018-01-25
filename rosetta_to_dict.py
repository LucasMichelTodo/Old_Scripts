#!/usr/bin/env python
# -*- coding: utf-8 -*-

rosetta = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Rosetta.txt", "r+") as file1:	
	for line in file1:
		rosetta[str(line.split("\t")[0].strip())] = {"old_names":line.split("\t")[1:], "old_fam":".", "annot":"."}
	for key, value in rosetta.iteritems():
		value["old_names"][-1] = value["old_names"][-1].strip()

with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Sorted_filtered.gff", "r") as file2:
	for line in file2:
		if line.split("\t")[8].split(";")[0].replace("ID=", "") in rosetta.keys():
			rosetta[line.split("\t")[8].split(";")[0].replace("ID=", "")]["annot"] = line.split("\t")[8].split(";")[1].strip().replace("description=", "")

with open("/home/lucas/ISGlobal/Chip_Seq/Transcripció_CSV/Taula_families_simplificada.csv", "r") as file3:
	for line in file3:
		for key, value in rosetta.iteritems():
			if line.split("\t")[0] in value["old_names"]:
				rosetta[key]["old_fam"] = line.split("\t")[1].strip()










