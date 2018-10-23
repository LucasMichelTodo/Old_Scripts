#!/usr/bin/env python

import sys
import re
import subprocess



subprocess.call("mkdir Tcd4_Epitope_Fastas", shell=True)

def parse(filein):
	i = 1
	with open(filein, "r+") as infile:
		for line in infile:
			epitope = line.strip().split()[0]
			ag = line.strip().split()[1]
			#organism = line.strip().split()[2]

			with open("Tcd4_Epitope_Fastas/"+re.sub('[\[\'\]]', '', epitope)+".fasta", "w+") as outfile:
				outfile.write(">"+"epitope_{}" .format(i) +" | "+ag+"\n")
			with open("Tcd4_Epitope_Fastas/"+re.sub('[\[\'\]]', '', epitope)+".fasta", "a+") as outfile:
				outfile.write(re.sub('[\[\'\]]', '', epitope))
			i += 1




if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		parse(element)
