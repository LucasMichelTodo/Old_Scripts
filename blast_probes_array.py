#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML



def blast_probes(probe_list):

    header = True
    hits = []

    with open(probe_list, "r+") as filein:
        # for line in filein:
        #     if header:
        #         header = False
        #         pass
        #     else:
        #         probe = line.split(",")[0].replace("/", ":")

        for line in filein:
            if line.startswith(">"):
                probe = line.strip().replace("/", ":").replace(">", "")

                cmd = ("blastn -query ./Missing_Probes_fastas/{}.fasta " +
                        "-subject /home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-37_Pfalciparum3D7_AnnotatedTranscripts.fasta " +
                        "-outfmt 5 -task blastn " +
                        "-dust no -soft_masking false " +
                        "> ./Missing_Probes_XMLs/blast_{}.xml") .format(probe,probe)

                print cmd
                subprocess.call(cmd, shell=True)

                result_handle = open("./Missing_Probes_XMLs/blast_{}.xml" .format(probe))
                blast_records = NCBIXML.parse(result_handle)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:

                            hits.append((alignment.title.split("|")[0].split()[1].strip(), hsp.bits))
                            #probe = i.replace("./blast_", "").replace(".fasta.xml", "")

                sorted_hits = sorted(hits, key=lambda x: x[1], reverse=True)

                if len(sorted_hits) < 1:
                    print probe
                    print "No hits!"
                else:
                    print probe

                    two_best_hits = []

                    for i in sorted_hits:
                        if (len(two_best_hits) >= 1) and (i[0].split(".")[0] in two_best_hits[0][0][0]):
                            if i[0] not in [x[0] for x in two_best_hits[0]]:
                                two_best_hits[0].append(i)

                        elif len(two_best_hits) < 2:
                            two_best_hits.append([i])

                        else:
                            pass

                    print two_best_hits

                print "---------------------------------------------------------------------"

                # with open("missing_probes_parsed_result_multitranscript.csv", "a+") as outfile:
                #     outfile.write(probe+"\t")
                #     for x in two_best_hits:
                #         outfile.write(x[0][0]+"\t"+str(x[0][1])+"\t")
                #     if len(two_best_hits[0]) > 1:
                #         first = True
                #         for i in two_best_hits[0]:
                #             if first:
                #                 outfile.write(str(i))
                #                 first = False
                #             else:
                #                 outfile.write(","+str(i))
                #
                #     outfile.write("\n")
                #
                # hits = []

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
