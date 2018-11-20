#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML

def parse_blastp(blast_xml):

    result_handle = open(blast_xml) #Run BLAST fom web interface first and download xml output (much faster)
    blast_records = NCBIXML.parse(result_handle)

    best_hits = []

    for record in blast_records:

        epitope_length = record.query.split(" | ")[1]
        #epitope_length = 9
        hits = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:

                hits.append((alignment.title, float(hsp.identities)/int(epitope_length)*100, hsp.expect, hsp.bits))
                #hits.append((alignment.title.split("|")[3], float(hsp.identities)/9*100, hsp.expect))
                sorted_hits = sorted(hits, key=lambda x: x[3], reverse=True)

        try:
            best_hits.append((record.query, sorted_hits[0]))
        except:
            best_hits.append((record.query, ""))
            #print "hopla!", record.query

    with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/Blasts/Microbiome/tcd4_Epitopes/tcd4_microbiome.csv", "r+") as infile:
        for line in infile:
            epitope = line.strip().split("\t")[0]
            # print line.strip().split("\t")
            # print len(line.strip().split("\t"))
            if len(line.strip().split("\t")) == 3:

                for i in best_hits:
                    if i[0].split(" | ")[0] == epitope:
                        blast_result = [str(x) for x in i[1]]
                        print line.strip()+"\t"+"\t"+"\t"+"\t"+"\t".join(blast_result)

            # print line.strip().split("\t")
            # print len(line.strip().split("\t"))
            else:
                for i in best_hits:
                    if i[0].split(" | ")[0] == epitope:
                        blast_result = [str(x) for x in i[1]]
                        print line.strip()+"\t"+"\t".join(blast_result)



if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        parse_blastp(i)
