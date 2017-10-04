#!/usr/bin/env python

import sys
import os
import subprocess
import itertools
import re

# Create function "run" to capture output and stderr of a "call".	

def run(cmd):
	call = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
	output, err = call.communicate()
	return (output, err)

# Provide a list of .bam files and make every possible 2 way comparison between them using macs2 bdgdiff. (macs2 callpeak is used previoulsy to comparison.)

def diff_peak_call(bams):
	if len(bams) < 2:
		sys.exit('Error! Less than 2 files provided.')

	subprocess.call("mkdir macs2_results_oldschool", shell=True)

	depths = {}
	names = set()

	for file in bams:
		names.add(file.split("_")[0]) # Create a set with the "names" of the alignments.

	for name in names:  #Callpeaks for each alignment vs it's input control.
		cmd = "macs2 callpeak -g 2.3e+7 -B -t {}_me_sort.bam -c {}_in_sort.bam -n {} --nomodel --extsize 117 --outdir ./macs2_results_oldschool" .format(name, name, name[2:])
		subprocess.call(cmd, shell=True)
		
		cmd = "egrep \"tags after filtering in treatment|tags after filtering in control\" ./macs2_results_oldschool/{}_peaks.xls" .format(name[2:]) #Get sequencing depth of the alignment.
		depth_info = run(cmd)[0]
		print depth_info

		num = re.compile("\d+")
		depths[name] = min(map(int,num.findall(depth_info)))
		print depths

	pairs = itertools.combinations(names, 2) # Create pairwise combinations
	for pair in pairs:
		# Run macs2 bdgdiff for each pairwise combination.
		cmd = "macs2 bdgdiff --t1 ./macs2_results_oldschool/{}_treat_pileup.bdg --c1 ./macs2_results_oldschool/{}_control_lambda.bdg --t2 ./macs2_results_oldschool/{}_treat_pileup.bdg --c2 ./macs2_results_oldschool/{}_control_lambda.bdg --d1 {} --d2 {} -g 150 -l 200 --o-prefix {} --outdir ./macs2_results_oldschool" \
		.format(pair[0], pair[0],
				 pair[1], pair[1], 
				 depths[pair[0]], depths[pair[1]], 
				 pair[0][2:]+"_"+pair[1][2:]) 
		
		subprocess.call(cmd, shell=True)



if __name__ == "__main__":
	diff_peak_call(sys.argv[1:])
