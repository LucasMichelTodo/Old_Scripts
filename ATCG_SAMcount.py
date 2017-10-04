#!/usr/bin/env python

# Import packages 
import sys
import os
from tqdm import tqdm
from Bio.Seq import Seq
import subprocess

# Fetch files from a given directory: hard coded

# fd_rawlist = os.listdir("/home/lucas/ISGlobal/TestSet/Prova/Results")
# fd_list = []
# for element in fd_rawlist:
# 	if element[-4:] == ".sam" and element[-6:-4] not in ("al", "un"):
# 		fd_list.append(element)

# # Set in/out paths
# inpath = "/home/lucas/ISGlobal/TestSet/Prova/Results/"

# COUNTING FUCTION
def run(cmd):
	call = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
	output, err = call.communicate()
	return (output, err)

def AAsamCount(samfile):
	mapping = {"-f4":"unmapped", "-F4":"mapped"}
	for i in mapping:
		mapsam = run("samtools view -h {} {}" .format(i,samfile))[0] # The output of "run" is a tupple with two elements: stdout and stderr
		with open(samfile.replace(".sam", "_{}.sam" .format(mapping[i])),"w+") as outfile:
			outfile.write(mapsam.strip("\n"))
		
		output = run("samtools view {}|cut -f10" .format(outfile.name))[0]
		seq = "" # Create empty sequence
		for line in tqdm(output): # Iterate over stdout of previous "run" (text file with only sequences from a sam file)
			seq += Seq(line.strip("\n")) # Concatenate sequences 
		

		with open(samfile.replace(".sam", "_AAcount.txt"),"a+") as result:
			result.write("COUNTED ATCG CONTENT FOR FILE: {}\n" .format(outfile.name)+"\n"													
						+"Total sequence length :"+str(len(seq))+"\n"
						+"\"A\" count: "+str(seq.count("A")/float(len(seq))*100)+"%"+"\n"
						+"\"T\" count: "+str(seq.count("T")/float(len(seq))*100)+"%"+"\n"
						+"\"G\" count: "+str(seq.count("G")/float(len(seq))*100)+"%"+"\n"
						+"\"C\" count: "+str(seq.count("C")/float(len(seq))*100)+"%"+"\n\n")


# Calling function for element in fd_list, from fetching from directory

# for element in fd_list:
# 	AAsamCount("{}{}" .format(inpath,element))

# Fetch files from console input

filenames = sys.argv[1:]
print filenames
for element in filenames:
	AAsamCount(element)






