#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML
from Bio import SeqIO


def blast_epitopes(epitope_fasta):

    header = False
    hits = []

    # with open(epitope_fasta, "r+") as filein:
    #     for line in filein:
    #         if line.startswith(">"):
    #             name = line.strip()[4:]
    #             ag = "epitopeB"
    #         else:
    #             epitope = line.strip()
    #             print name, ag, epitope, str(len(epitope))

    record = SeqIO.read(epitope_fasta, "fasta")
    name = str(record.id[3:])
    ag = "epitopeB"
    epitope = str(record.seq)
    #print name, ag, epitope, str(len(epitope))

    cmd = ("blastp -query {} " +
            "-db /home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Blasts/microbiome_fastas.fasta " +
            "-outfmt 5 -task blastp " +
            "-matrix PAM30 " +
            "-num_threads 4 -evalue 1000 "
            #"-dust no -soft_masking false " +
            "> /media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Exposed/B_predictions/Blasts/blast_{}.xml") .format(epitope_fasta.replace("./", ""),name)

    #print cmd
    subprocess.call(cmd, shell=True)

    result_handle = open("/media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Exposed/B_predictions/Blasts/blast_{}.xml" .format(name))

    try:
        blast_records = NCBIXML.parse(result_handle)

        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:

                    hits.append((alignment.title.split("|")[5], hsp.bits, float(hsp.identities)/len(epitope)*100, hsp.expect))
                    #probe = i.replace("./blast_", "").replace(".fasta.xml", "")

        sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)

        if len(sorted_hits) < 1:
            print "No hits!"
        else:

            pass
        #     for i in sorted_hits[0:5]:
        #         print i
        #         print "\n"
        #
        # print "---------------------------------------------------------------------"

        with open("/media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Exposed/B_predictions/blasted_microbiome.csv", "a+") as outfile:
            outfile.write(name+"\t"+epitope+"\t"+ag+"\t")
            if len(sorted_hits) > 0:
                outfile.write(str(sorted_hits[0][0])+"\t"+str(sorted_hits[0][2])+"\t"+str(sorted_hits[0][3]))

            outfile.write("\n")

    except:
        print "no results for: ", name


if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        blast_epitopes(i)
