#!/usr/bin/env python

dic1 = { "mRNA":"SO:0000234",
		"gene":"SO:0000704",
		"exon":"SO:0000147",
		"CDS":"SO:0000316",
		"ncRNA":"SO:0000655",
		"rRNA":"SO:0000252",
		"snoRNA":"SO:0000275",
		"snRNA":"SO:0000274",
		"three_prime_UTR":"SO:0000205",
		"tRNA":"SO:0000253"}

dic2 = { "mRNA":"protein_coding",
		"gene":"protein_coding",
		"exon":"protein_coding",
		"CDS":"protein_coding",
		"ncRNA":"non_coding",
		"rRNA":"rRNA",
		"snoRNA":"non_coding",
		"snRNA":"non_coding",
		"three_prime_UTR":"non_coding",
		"tRNA":"Mt_tRNA"}

with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", "r+") as file1:
		for line in file1:
			if line.startswith("#"):
				print line.strip()
			else:
				line_list = line.split("\t")
				#line_list[2] = dic1[line.split()[2]]
				print "\t".join(line_list).strip()+";biotype="+dic2[line.split()[2]]
