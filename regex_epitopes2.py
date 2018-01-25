#!/usr/bin/env python

import os
import sys
import re
from Bio import SeqIO
import pybedtools as py

epitopes_t = []
epitopes_b = []

gff = py.BedTool("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like.gff")

# Create a list with epitopes from file and convert them to regex. Change this file to change epitopes to search.
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/tcd8_epitopes_to_cross.txt", "r+") as file:
	for line in file:
		regex = re.compile(re.escape(str(line.strip()))) 
		epitopes_t.append(regex)

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/bcell_epitopes_to_cross.txt", "r+") as file:
	for line in file:
		regex = re.compile(re.escape(str(line.strip()))) 
		epitopes_b.append(regex)

#Create header for results (ensures file is written from scratch if re-run).
with open("../../Results_tcell.txt", "w+") as file: 
	file.write("cluster"+"\t"+"sequence"+"\t"+"epitope"+"\t"+"number of matches"+"\t"+"annotation"+"\n")

with open("../../Results_bcell.txt", "w+") as file: 
	file.write("cluster"+"\t"+"sequence"+"\t"+"epitope"+"\t"+"number of matches"+"\t"+"annotation"+"\n")


def check_epitopes(fasta_file):

	if os.stat(fasta_file).st_size == 0: # Break funtion if fasta file is empty.
		return

	seq = ""
	with open(fasta_file, "r+") as file: # fetch sequence.
		for line in file:
			if not line.startswith(">"):
				seq += str(line)
			else:
				prot = line.strip().replace(">","").split("/")[0]



	for epitope in epitopes_t:
		result = re.findall(epitope, seq) # find matches in sequence.

		if result: 
			
			for record in gff:
				if record[8].split(";")[0].replace("ID=","") == prot:
					anot = record[8]

			if 'anot' not in locals():
				anot = "."

			print anot

			print fasta_file+"\n"
			print prot+"\n" 
			print result[0]+"("+str(len(result))+")\n"
			print "\n"+"-----------------------------------------------------------------"+"\n"
			with open("../../Results_tcell.txt", "a+") as file:
				file.write(fasta_file.replace("_consensus.fasta", "")+"\t"+prot+"\t"+result[0]+"\t"+str(len(result))+"\t"+anot+"\n") # Write results.

	for epitope in epitopes_b:
			result = re.findall(epitope, seq) # find matches in sequence.

			if result:
				
				for record in gff:
					if record[8].split(";")[0].replace("ID=","") == prot:
						anot = record[8]

				if 'anot' not in locals():
					anot = "."

				print anot

				print fasta_file+"\n"
				print prot+"\n" 
				print result[0]+"("+str(len(result))+")\n"
				print "\n"+"-----------------------------------------------------------------"+"\n"
				with open("../../Results_bcell.txt", "a+") as file:
					file.write(fasta_file.replace("_consensus.fasta", "")+"\t"+prot+"\t"+result[0]+"\t"+str(len(result))+"\t"+anot+"\n") # Write results.

if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:	
		check_epitopes(file)
