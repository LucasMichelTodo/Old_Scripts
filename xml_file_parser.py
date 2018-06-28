#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO

Entrez.email = "lucas.michel@isglobal.org"     # Always tell NCBI who you are


print "Antigen\tHit_1\tScore_1\tHit_2\tScore_2"
#Parse blast_xml output
def parse_blastp(blast_xml):

    result_handle = open(blast_xml) #Run BLAST fom web interface first and download xml output (much faster)
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        hits = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:

                hits.append((":".join(alignment.title.strip().split("|")[3:5]).replace(">gi", "").strip(), hsp.bits))

                sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)


        print_list = [record.query]
        for x in sorted_hits[0:2]:
            print_list.append(str(x[0]))
            print_list.append(str(x[1]))

        print "\t".join(print_list)

    result_handle.close()

if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        parse_blastp(file)
