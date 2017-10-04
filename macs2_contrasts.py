#!/usr/bin/env python

import sys
import os
import subprocess
import itertools
import argparse

parser = argparse.ArgumentParser(description='Peak calling for contrasts')

parser.add_argument('bams', type=str, nargs='+', help='Bam files to be called.')
parser.add_argument('-f','--fe', action="store", type=float, default = 1.0, dest="fe", help='Fold enrichment cutoff (float).')
parser.add_argument('-q','--qval', action="store", type=float, default = 0.05, dest="qval", help='q-value cutoff (float).')
parser.add_argument('-b', '--broad', action="store_true", default = False, dest="broad", help='Call peaks using "mac2 calpeak" --broad option.')
parser.add_argument('-n', '--dirname', action="store", type=str, default = "macs2_constrast", dest="dirname", help='Name of the folder that will be created to hold the results.')

# Create function "run" to capture output and stderr of a "call".	

def run(cmd):
	call = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
	output, err = call.communicate()
	return (output, err)

# Provide a list of .bam files and make every possible 2 way comparison between them using macs2 bdgdiff. (macs2 callpeak is used previoulsy to comparison.)

def diff_peak_call(bams, fe, qval, broad):
	if len(bams) < 2:
		sys.exit('Error! Less than 2 files provided.')

	subprocess.call("mkdir {}" .format(args.dirname), shell=True)
	
	if broad: # Wether tu use "broad" mode on macs2
		usebroad = "--broad"
	else:
		usebroad = ""
	
	names = set()

	for file in bams:
		names.add(file.split("_")[0]) # Create a set with the "names" of the alignments.

	pairs = itertools.combinations(names, 2)
	for pair in pairs:  #Callpeaks for each pair of alignments (each contrast).
		cmd = "macs2 callpeak -g 2.3e+7 -f BAMPE -t {}_me_sort.bam -c {}_me_sort.bam -n {} --nomodel --extsize 117 --fe-cutoff {} -q {} {} --outdir ./{}" \
		.format(pair[0], pair[1], pair[0][2:]+"_"+pair[1][2:], fe, qval, usebroad, args.dirname)
		subprocess.call(cmd, shell=True)
		
		cmd = "macs2 callpeak -g 2.3e+7 -f BAMPE -t {}_me_sort.bam -c {}_me_sort.bam -n {} --nomodel --extsize 117 --fe-cutoff {} -q {} {} --outdir ./{}" \
		.format(pair[1], pair[0], pair[1][2:]+"_"+pair[0][2:], fe, qval, usebroad, args.dirname)
		subprocess.call(cmd, shell=True)
		

if __name__ == "__main__":
	args = parser.parse_args()
	diff_peak_call(args.bams, args.fe, args.qval, args.broad)
