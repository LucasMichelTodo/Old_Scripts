#!/usr/bin/env python

import sys
import os
import subprocess

#Create dirs
subprocess.call("mkdir ./cluster_withref", shell=True)
subprocess.call("mkdir ./cluster_noref", shell=True)
subprocess.call("mkdir ./cluster_withref/fastas", shell=True)
subprocess.call("mkdir ./cluster_withref/msa", shell=True)
subprocess.call("mkdir ./cluster_noref/fastas", shell=True)
subprocess.call("mkdir ./cluster_noref/msa", shell=True)

# Call all perl scripts sequentially
subprocess.call("perl parse_tcruzi-clstr_01.pl", shell=True)
subprocess.call("perl muscle_tcruzi-clustr.pl", shell=True)
subprocess.call("perl sort_msa.pl", shell=True)

# Erase a file called "_sorted.aln" that is accidentally created.
try:
    subprocess.call("rm cluster_withref/msa/_sorted.aln", shell = True) #Sha darreglar!!!!
    print "Erasing empty files!"
except:
    print "Nothing to erase!"

# Call all python scripts sequentially
subprocess.call("entropy.py ./cluster_withref/msa/*\_sorted.aln", shell=True)
subprocess.call("regex_epitopes2.py ./cluster_withref/msa/*\_consensus.fasta", shell=True)
