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

subprocess.call("perl parse_tcruzi-clstr_01.pl", shell=True)
subprocess.call("perl muscle_tcruzi-clustr.pl", shell=True)
subprocess.call("perl sort_msa.pl", shell=True)
#subprocess.call("rm ./cluster_withref/msa/_sorted.aln", shell=True)

subprocess.call("rm ./cluster_withref/msa/_sorted.aln ./cluster_withref/msa/_sorted.msa") #Sha darreglar!!!!
subprocess.call("entropy.py ./cluster_withref/msa/*\_sorted.aln", shell=True)

