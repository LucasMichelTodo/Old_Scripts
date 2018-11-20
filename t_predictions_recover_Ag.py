#!/usr/bin/env python

import sys
reference = {}
annotation = {}

## Put your reference file here:
reference_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_consensus00.fasta"
annotation_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_proteomes.fasta"

def recover_origin_ag(epitopes_file):

    with open(reference_file, "r+") as ref:
        for line in ref:
            if line.startswith(">"):
                prot = line.strip().split("/")[0].replace(">", "")
            else:
                seq = line.strip()
                reference[prot] = seq

    with open(annotation_file, "r+") as ano_ref:
        for line in ano_ref:
            if line.startswith(">"):
                annotation[line.strip().split(" | ")[0].replace(">", "")] = [x.replace("gene_product=", "") for x in line.strip().split(" | ") if x.startswith("gene_product=")][0]

    # for key, value in annotation.iteritems():
    #     print key
    #     print value

    with open(epitopes_file, "r+") as filein:
        for line in filein:
            epitope = line.strip()
            for ag, seq in reference.iteritems():
                if epitope in seq:
                    if ag in annotation.keys():
                        print "{}\t{}\t{}" .format(epitope, ag, annotation[ag])
                    else:
                        print "{}\t{}\t." .format(epitope, ag)


if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		recover_origin_ag(file)
