#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO

#Parse blast_xml output
def parse_blastp(blast_xml):
    x = 0
    i = 0
    result_handle = open(blast_xml) #Run BLAST fom web interface first and download xml output (much faster)
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        hits = []

        for alignment in record.alignments:
            for hsp in alignment.hsps:

                prot = record.query_id
                seq = record.query
                e_val = hsp.expect
                score = hsp.bits
                identities = hsp.identities
                identity = round(float(hsp.identities)/hsp.align_length*100, 2)

                hits.append((prot, seq, alignment.title, e_val, score, identities, identity))

        sorted_hits = sorted(hits, key=lambda x: x[4], reverse=True)

        try:
            if  sorted_hits[0][3] < 10 and sorted_hits[0][6] > 70:
                print "Query:   ", sorted_hits[0][0]
                print "Query Annot:   ", sorted_hits[0][1]
                print "Hit:   ", sorted_hits[0][2]
                print "E-val:   ", sorted_hits[0][3]
                print "Score:   ", sorted_hits[0][4]
                print "Identities:   ", sorted_hits[0][5]
                print "Identity:   ", sorted_hits[0][6],"%"
                print "----------------------------------"
        except:
            print "No hits"
            print "----------------------------------"

    result_handle.close()




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        parse_blastp(file)
