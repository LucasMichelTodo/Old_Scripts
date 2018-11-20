#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO

#Parse blast_xml output
def parse_blastp(blast_xml):
    x = 1
    result_handle = open(blast_xml) #Run BLAST fom web interface first and download xml output (much faster)
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        print x
        x += 1
        i = 1
        if record.alignments:
            for alignment in record.alignments:


                for hsp in alignment.hsps:

                    print "******************Alignment {}.{} **************************" .format(x,i)
                    print record.query_id
                    print record.query
                    print hsp.expect
                    print hsp.bits
                    print hsp.identities
                    print round(float(hsp.identities)/hsp.align_length*100, 2)

                    i += 1
        else:
            print record.query, "no hits"

    result_handle.close()


if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        parse_blastp(file)
