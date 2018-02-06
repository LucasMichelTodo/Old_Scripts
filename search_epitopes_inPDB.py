#!/usr/bin/env python

import regex
import sys
from tqdm import tqdm

epitopes = []
epitopes_fromfile = []
fasta = {}

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/tcd8_epitopes_to_cross.txt", "r+") as file:
	for line in file:
		epitopes.append(regex.compile("(?b)"+"("+line.strip()+")"+"{e<=1}"))

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/bcell_epitopes_to_cross.txt", "r+") as file:
	for line in file:
		epitopes.append(regex.compile("(?b)"+"("+line.strip()+")"+"{e<=1}"))

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/PDB_fastas/pdb.fasta", "r+") as file1:
    for line in file1:
        if line.startswith(">"):
            current_gene = line.strip()
            fasta[current_gene] = ""
        else:
            fasta[current_gene] += line.strip()

hits = []
for epi in tqdm(epitopes):
    for key, value in fasta.items():
        if epi.search(value):
            hits.append(key[1:5])

for i in set(hits):
    print i

def search_epitopes(epitopes_file):
    with open(epitopes_file, "r+") as file:
        firstline = True
    	for line in file:
            if firstline:
                firstline = False
            else:
                epitopes_fromfile.append(regex.compile("(?b)"+"("+line.split()[2]+")"+"{e<=1}"))

    for epi in epitopes_fromfile:
        for key, value in fasta.items():
            if epi.search(value):
                print key
                print epi.pattern.replace("(?b)", "").replace("{e<=1}", "")
                print epi.search(value)
                print "----------------------------------------------------------------------------------------------"


if __name__ == "__main__":
    filenames= sys.argv[1:]
    for file in filenames:
        search_epitopes(file)
