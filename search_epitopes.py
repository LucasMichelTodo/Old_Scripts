#!/usr/bin/env python

import os
import sys
import re
from Bio import SeqIO
import pybedtools as py

def fastatodict(fasta_file):

    fasta = {}
    seq = ""
    prot = ""

    with open(fasta_file, "r+") as file1:
        for line in file1:
            if line.startswith(">"):
                fasta[prot] = seq
                seq = ""
                prot = line.strip()
            else:
                seq += line.strip()

                fasta[prot] = seq

    return fasta

def search_epitopes(epitopes_file, fasta_file):

	epitopes = []

	gff = py.BedTool("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like.gff")

	# Create a list with epitopes from file and convert them to regex. Change this file to change epitopes to search.
	with open(epitopes_file, "r+") as file:
		for line in file:
			regex = re.compile(re.escape(str(line.strip())))
			if len(regex.pattern) > 1:
				epitopes.append(regex)


	#Create header for results (ensures file is written from scratch if re-run).
	results_file = "search_epitopes_"+epitopes_file.split(".")[0]+"_"+fasta_file.split(".")[0]+".txt"

	with open(results_file, "w+") as file:
		file.write("cluster"+"\t"+"sequence"+"\t"+"epitope"+"\t"+"number of matches"+"\t"+"annotation"+"\n")

	if os.stat(fasta_file).st_size == 0: # Break funtion if fasta file is empty.
		return

	fasta = fastatodict(fasta_file)

	for prot, seq in fasta.iteritems():
		if len(prot) > 0:
			for epitope in epitopes:
				result = re.findall(epitope, seq) # find matches in sequence.


				if len(result) > 0:

					for record in gff:
						if record[8].split(";")[0].replace("ID=","") == prot:
							anot = record[8]

					if 'anot' not in locals():
						anot = "."

					with open(results_file, "a+") as file:
						file.write("unknown_cluster"+"\t"+prot+"\t"+result[0]+"\t"+str(len(result))+"\t"+anot+"\n") # Write results.


if __name__ == "__main__":
	filenames= sys.argv[1:]
	search_epitopes(filenames[0], filenames[1])
