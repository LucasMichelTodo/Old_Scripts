#!/usr/bin/env python

gff = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Sorted_filtered.gff", "r") as file1:
	for line in file1:
		gff[line.plit()[8].split(";")[1].replace("ID=", "")]={"start":line.plit()[3], "stop":line.plit()[4], "chrom":line.plit()[0], "strain":line.plit()[6], "anot":line.plit()[8]}

print gff[1:3]