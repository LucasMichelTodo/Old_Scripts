#!/usr/bin/env python

import sys

probe_dict = {}

with open("/home/lucas/ISGlobal/Arrays/Array_Annotation/probes.fasta", "r+") as file_ref:
    for line in file_ref:
        if line.startswith(">"):
            probe = line.strip().replace(">", "")
        else:
            seq = line.strip()
            probe_dict[probe] = seq

probes_noseq = []

with open("/home/lucas/ISGlobal/Arrays/Array_Annotation/missing_probes.txt", "r+") as filein:
    for line in filein:
        probe = line.strip()
        if probe in probe_dict.keys():
            print ">"+probe
            print probe_dict[probe]
            with open("/home/lucas/ISGlobal/Arrays/Missing_Probes_fastas/{}.fasta" .format(probe.replace("/", ":")), "w+") as outfile:
                outfile.write(">"+probe+"\n"+probe_dict[probe])
        else:
            probes_noseq.append(probe)

# for i in probes_noseq:
#     print i
