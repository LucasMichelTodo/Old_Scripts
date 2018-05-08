#!/usr/bin/env python

fasta = {}
seq = ""

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/ncbi_pdb_tcruzi.fasta", "r+") as file1:
    for line in file1:
        if line.startswith(">"):
            line_split = line.strip().split("|")
            fasta[line_split[1]] = {"gi":line_split[1], "pdb":line_split[3], "chain":line_split[4][0], "annot":line_split[4][11:], "seq":""}
            current_gi = line_split[1]
        elif len(line.strip()) < 1:
            print current_gi
            print seq
            fasta[current_gi]["seq"] = seq
            seq = ""
        else:
            seq += line.strip()

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/ncbi_pdb_tcruzi.txt", "w+") as file2:
    for key, value in fasta.iteritems():
        file2.write(value["gi"]+"\t"+value["pdb"]+"\t"+value["chain"]+"\t"+value["annot"]+"\t"+value["seq"]+"\n")
