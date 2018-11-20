#!/usr/bin/env python

import sys
import subprocess
from blast_epitopes2 import blast_epitopes
from Bio.Blast import NCBIXML

## Pass a plain txt file with one epitope per line as an input.

epitopes = []
def from_list_to_fasta(epitopes_list):

    subprocess.call("mkdir Fastas", shell=True)
    i = 1

    with open(epitopes_list, "r+") as epitopeslist:
        for line in epitopeslist:
            epitope = line.strip()
            epitopes.append(epitope)
            with open("./Fastas/"+epitope+".fasta", "w+") as fastafile:
                fastafile.write(">Epitope_{} | {}" .format(i, len(epitope))+"\n"+epitope)
                i+=1

    for i in epitopes:
        blast_epitopes("./Fastas/"+i+".fasta")




if __name__ == "__main__":
    filename= sys.argv[1]
    from_list_to_fasta(filename)
