#!/usr/bin/env python

import sys
import re

epitopes = []

# Create a list with epitopes from file and convert them to regex.
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/tcd8epitopes_to_cross.txt", "r+") as file:
	for line in file:
		regex = re.compile(re.escape(str(line.strip()))) 
		epitopes.append(regex)

# Create header for results (ensures file is written from scratch if re-run).
with open("../../Results.txt", "w+") as file: 
	file.write("cluster"+"\t"+"sequence"+"\t"+"epitope"+"\t"+"number of matches"+"\n")

def check_epitopes(fasta_file):

	with open(fasta_file, "r+") as file: # fetch sequence.
		for line in file:
			if not line.startswith(">"):
				seq = str(line)
			else:
				prot = line.strip()

	for epitope in epitopes:
		result = re.findall(epitope, seq) # find matches in sequence.
		if result:
			print fasta_file+"\n"
			print prot+"\n" 
			print result[0]+"("+str(len(result))+")\n"
			print "\n"+"-----------------------------------------------------------------"+"\n"
			with open("../../Results.txt", "a+") as file:
				file.write(fasta_file.replace("_consensus.fasta", "")+"\t"+prot+"\t"+result[0]+"\t"+str(len(result))+"\n") # Write results.

if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:	
		check_epitopes(file)
