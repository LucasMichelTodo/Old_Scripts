#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML



def blast_probes(probe_list):

    # with open("Result.csv", "w+") as outfile:
    #     outfile.write("\t".join(["Probe", "Hit", "Score"]))
    #     outfile.write("\n")
    #
    for i in probe_list:
        hits = []
    #     print i
    #     cmd = "blastn -query {} -subject /home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-30_Pfalciparum3D7_AnnotatedTranscripts.fasta -outfmt 5 -task blastn -max_target_seqs 2 -max_hsps 2 > blast_{}.xml" .format(i,i.replace("./", ""))
    #     subprocess.call(cmd, shell=True)

        result_handle = open(i)
        #result_handle = open("blast_{}.xml" .format(i.replace("./", "")))
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:

                    hits.append((alignment.title.split("|")[0].split()[1].strip(), hsp.bits))
                    probe = i.replace("./blast_", "").replace(".fasta.xml", "")

        sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)

        if len(sorted_hits) < 1:
            probe = i.replace("./blast_", "").replace(".fasta.xml", "")
            print probe
            print "No hits!"

        # print probe
        # print sorted_hits[0:2]
        #
        # with open("parsed_result_bits2.csv", "a+") as outfile:
        #     outfile.write(probe+"\t")
        #     for x in sorted_hits[0:2]:
        #         outfile.write(x[0]+"\t"+str(x[1])+"\t")
        #     outfile.write("\n")





if __name__ == "__main__":
    filenames = sys.argv[1:]
    blast_probes(filenames)
