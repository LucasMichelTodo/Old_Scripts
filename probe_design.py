#!/usr/bin/env python

import sys
from Bio import AlignIO


def  search_probe(align_file):
        alignment = AlignIO.read(align_file, "clustal")
        probe_list = []

        print str(align_file)

        for i in range(len(alignment._records[0])-60):
            probe = alignment._records[0].seq[i:i+60]
            probe_score = alignment._star_info[i:i+60].count("*")
            if "-" not in probe._data:
                probe_list.append((probe, probe_score))
        
            #print probe
            #print probe_score

        sorted_probes = sorted(probe_list, key=lambda x: x[1])
        for i in sorted_probes[1:20]:
            print i[0]._data, i[1]

        print "\n **************************************************************************************************** \n"




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        search_probe(i)
