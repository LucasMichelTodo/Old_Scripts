#!/usr/bin/env python

import re

## Read reference fasta file and epitopes file
fasta_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/data/proteomes/TriTrypDB-38_TcongolenseIL3000_AnnotatedProteins.fasta"
epitopes_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/IEDB_recommended/Tcd4_predictions/tcd4_cores_and_consensus.csv"
fasta = {}


## Parse reference fasta and convert it into a dictionary
with open(fasta_file, "r+") as infile1:
    for line in infile1:
        if line.startswith(">"):
            prot = line.strip()
            fasta[prot] = ""
        else:
            fasta[prot] += line.strip()


## Parse epitopes file and compile a regex for each epitope and put it in a list
epitopes = []
with open(epitopes_file, "r+") as infile:

    for line in infile:
        regex = re.compile(line.split()[1].strip())
        epitopes.append(regex)


## For each epitope(regex) search it into the fasta dictionary and create a dict to store fasta entries for each epitope:
result_dict = {}
for i in epitopes:
    for key, item in fasta.iteritems():
        result = re.findall(i, item)
        if result:
            ano = [x.replace("gene_product=", "") for x in key.split(" | ") if x.startswith("gene_product=")][0]
            result_dict.setdefault(i.pattern,[]).append(ano)
            # print i.pattern
            # print key
            # print item

## Print Results
for ep in epitopes:
    try:
        print ep.pattern+"\t"+result_dict[ep.pattern][0]
    except:
        print ep.pattern
