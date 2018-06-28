#!/usr/bin/env python

rosetta = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", "r+") as file1:
	for line in file1:
		rosetta[str(line.split("\t")[0].strip())] = {"old_refs":line.split("\t")[1:]}
	for key, value in rosetta.iteritems():
		value["old_refs"][-1] = value["old_refs"][-1].strip()

with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", "r+") as file2:
	for line in file2:
		if line.startswith("#"):
			pass
		elif line.split()[2] == "gene":
			line_split = line.strip().split("\t")[8].split(";")
			if  line_split[0].replace("ID=", "") in rosetta.keys():
				rosetta[line_split[0].replace("ID=", "")]["annot"] = line_split[1].replace("description=", "")
		else:
			pass

with open("/home/lucas/ISGlobal/Gen_Referencies/gene_names.txt", "r+") as file3:
	header = True
	for line in file3:
		if header:
			pass
			header = False
		else:
			if line.strip().split()[0] in rosetta.keys():
				if line.strip().split("\t")[4] == "N/A":
					rosetta[line.strip().split("\t")[0]]["name"] = "NA"
				else:
					rosetta[line.strip().split("\t")[0]]["name"] = line.strip().split("\t")[4]

print rosetta
