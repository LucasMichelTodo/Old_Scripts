#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO

Entrez.email = "lucas.michel@isglobal.org"     # Always tell NCBI who you are

#Parse blast_xml output
def parse_blastp(blast_xml):
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


                hits.append((prot, seq, alignment.title.split(">")[0], e_val, score, identities, identity))

                # if float(hsp.identities)/float(hsp.align_length) > 0.80:
                #
                #     #Search Entrez query sequence to check if length of the hit alignment covers 90% of it:
                #     handle = Entrez.efetch(db="nucleotide", id=record.query_id, rettype="fasta", retmode="text")
                #     fasta = handle.read()
                #     handle.close()
                #     prot = fasta.split("\n")[0]
                #     seq = "".join(fasta.strip().split("\n")[1:])
                #
                #     if float(hsp.align_length)/float(len(seq)) > 0.70:

                # print('****Alignment {}****') .format(i)
                # print prot
                # print seq
                # print('sequence:', alignment.title.split(">")[0])
                # print('length:', alignment.length)
                # print('e value:', hsp.expect)
                # print('score:', hsp.score)
                # print('bitscore:', hsp.bits)
                # print hsp.identities
                # print hsp.align_length
                # print "-------------------"
                # print(hsp.query)
                # print "-------------------"
                # print(hsp.sbjct)
                # print "-------------------"
                # print hsp.frame
                # print hsp.query_start
                # i += 1
                        # with open("fasta_withPDBmodel.fasta", "a+") as outfile:
                        #     outfile.write(prot+"\n"+seq+"\n")

        sorted_hits = sorted(hits, key=lambda x: x[4], reverse=True)

        print "Hits: ",len(sorted_hits)

        strains = []
        for i in sorted_hits:
            organism = i[2].split(" | ")
            organism = [x for x in organism if x.startswith('organism=')][0]

            if organism not in strains:
                # print "Query:   ", i[0]
                # print "Query Annot:   ", i[1]
                print ">"+i[2]
                # print "E-val:   ", i[3]
                # print "Score:   ", i[4]
                # print "Identities:   ", i[5]
                # print "Identity:   ", i[6],"%"
                # print "----------------------------------"
                strains.append(organism)


    result_handle.close()




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        parse_blastp(file)
