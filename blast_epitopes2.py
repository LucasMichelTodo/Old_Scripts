#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML



def blast_epitopes(epitope_fasta):

    header = False
    hits = []

    with open(epitope_fasta, "r+") as filein:
        for line in filein:
            if line.startswith(">"):
                name = line.split(" | ")[0].strip().replace(">", "")
                ag = "epitopeB"
            else:
                epitope = line.strip()

    #print name, ag, epitope, str(len(epitope))+"\n"


    cmd = ("blastp -query {} " +
            "-db /home/lucas/ISGlobal/Brucei/aat_vaccine/data/Microbiome/all_microbiome.fasta " +
            "-outfmt 5 -task blastp " +
            "-matrix PAM30 " +
            "-num_threads 4 -evalue 10000000 -seg no "
            #"-dust no -soft_masking false " +
            "> /home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/Blasts/Microbiome/blast_{}.xml") .format(epitope_fasta.replace("./", ""),name)

    #print cmd
    subprocess.call(cmd, shell=True)

    result_handle = open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/Blasts/Microbiome/blast_{}.xml" .format(name))

    try:
        blast_records = NCBIXML.parse(result_handle)

        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:

                    hits.append((alignment.title, hsp.bits, float(hsp.identities)/len(epitope)*100, hsp.expect))
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

        with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/Blasts/Microbiome/tcd4_microbiome.csv", "a+") as outfile:
            outfile.write(name+"\t"+epitope+"\t"+ag+"\t")
            if len(sorted_hits) > 0:
                outfile.write(str(sorted_hits[0][0])+"\t"+str(sorted_hits[0][2])+"\t"+str(sorted_hits[0][3]))

            outfile.write("\n")

    except:
        print "no results for: ", name



    #result_handle = open(i)
    #result_handle = open("blast_{}.xml" .format(i.replace("./", "")))
    # blast_records = NCBIXML.parse(result_handle)
    # x = 0
    # for record in blast_records:
    #     for alignment in record.alignments:
    #         while x < 1:
    #             print('****Alignment {}****') .format(i)
    #             for hsp in alignment.hsps:
    #                 print ('query:', record.query_id)
    #                 print('sequence:', alignment.title)
    #                 print('length:', alignment.length)
    #                 print('e value:', hsp.expect)
    #                 print hsp.identities
    #                 print hsp.align_length
    #                 print "Score: {}" . format(hsp.score)
    #                 print "-------------------"
    #                 print(hsp.query)
    #                 print "-------------------"
    #                 print(hsp.sbjct)
    #                 print "-------------------"
    #                 print hsp.frame
    #                 print hsp.query_start
    #
    #                 # probe = i.replace("./blast_", "").replace(".fasta.xml", "")
    #                 # print probe
    #                 x += 1




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        blast_epitopes(i)
