#!/usr/bin/env python

import sys
reference = {}

def recover_origin_ag(epitopes_file):

    with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_27_02_18/fastas/all_05_consensus.fasta", "r+") as ref:
        for line in ref:
            if line.startswith(">"):
                prot = line.strip().split("/")[0].replace(">", "")
            else:
                seq = line.strip()
                reference[prot] = seq

    with open(epitopes_file, "r+") as filein:
        for line in filein:
            epitope = line.strip()
            for ag, seq in reference.iteritems():
                if epitope in seq:
                    print "{}\t{}" .format(epitope, ag)


if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		recover_origin_ag(file)
