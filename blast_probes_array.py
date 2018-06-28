#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML



def blast_probes(probe_list):

    header = True
    hits = []

    with open(probe_list, "r+") as filein:
        for line in filein:
            if header:
                header = False
                pass
            else:
                probe = line.split(",")[0].replace("/", ":")

                cmd = ("blastn -query ./Probes_fastas/{}.fasta " +
                        "-subject /home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-30_Pfalciparum3D7_AnnotatedTranscripts.fasta " +
                        "-outfmt 5 -task blastn -max_target_seqs 2 " +
                        "-max_hsps 2 > blast_{}.xml") .format(probe,probe)

                print cmd
                subprocess.call(cmd, shell=True)

                result_handle = open("blast_{}.xml" .format(probe))
                blast_records = NCBIXML.parse(result_handle)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:

                            hits.append((alignment.title.split("|")[0].split()[1].strip(), hsp.bits))
                            probe = i.replace("./blast_", "").replace(".fasta.xml", "")

                sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)

                if len(sorted_hits) < 1:
                    print probe
                    print "No hits!"
                else:
                    print probe
                    print sorted_hits[0:2]

                print "---------------------------------------------------------------------"

        # result_handle = open(i)
        # #result_handle = open("blast_{}.xml" .format(i.replace("./", "")))
        # blast_records = NCBIXML.parse(result_handle)
        # for record in blast_records:
        #     for alignment in record.alignments:
        #         print('****Alignment {}****') .format(i)
        #         for hsp in alignment.hsps:
        #             print ('query:', record.query_id)
        #             print('sequence:', alignment.title)
        #             print('length:', alignment.length)
        #             print('e value:', hsp.expect)
        #             print hsp.identities
        #             print hsp.align_length
        #             print "Score: {}" . format(hsp.score)
        #             print "-------------------"
        #             print(hsp.query)
        #             print "-------------------"
        #             print(hsp.sbjct)
        #             print "-------------------"
        #             print hsp.frame
        #             print hsp.query_start
        #
        #             probe = i.replace("./blast_", "").replace(".fasta.xml", "")
        #             print probe




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        blast_probes(i)
