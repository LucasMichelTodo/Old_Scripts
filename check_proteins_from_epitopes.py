#!/usr/bin/env python

import sys
import subprocess

fasta = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/TriTrypDB-35_TcruziCLBrenerEsmeraldo-like_AnnotatedProteins.fasta", "r+") as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:].split(":")[0]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = ""
            continue
        sequence = line
        fasta[active_sequence_name] += sequence

def check_prot(result_file):
	new_fasta = {}
	with open(result_file, "r+") as file:
		title = True
		for line in file:
			if title:
				title = False
			else:
				prot = line.split("\t")[1]
				seq = fasta[prot]
				if prot not in new_fasta:
					new_fasta[prot] = seq
	
	with open(result_file.replace(".txt", ".fasta"), "a+") as out_file:
		for key in new_fasta:
			out_file.write(">"+key+"\n"+new_fasta[key]+"\n")

	cmd = "signalp "+result_file.replace(".txt", ".fasta")+" > "+result_file.replace(".txt", "_signalp.txt")
	subprocess.call(cmd, shell=True)

	cmd = "targetp "+result_file.replace(".txt", ".fasta")+" > "+result_file.replace(".txt", "targetp.txt")
	subprocess.call(cmd, shell=True)

	#-------------------------------------------------------------------------------------------------------------------------------------------------------

	signal_or_not = {}
	with open(result_file.replace(".txt", "_signalp.txt"), "r+") as file_two:
		for line in file_two:
			if line.startswith("#"):
				pass
			else:
				signal_or_not[line.split()[0]] = line.split()[9]

	with open(result_file, "r+") as file:
		title = True
		for line in file:
			if title:
				title = False
			else:
				prot = line.split("\t")[1]
				with open(result_file.replace(".txt", "_tableSP.txt"), "a+") as outfile_two:
					outfile_two.write(line.strip()+"\t"+signal_or_not[prot]+"\n")





if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:	
		check_prot(file)






