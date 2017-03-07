#!/usr/bin/env python

from tqdm import tqdm
from Bio.Seq import Seq

def get_ACGTcontent_from_FASTA(filename):
	fd= open(filename)
	seq=""
	for line in tqdm(fd): 
		if not line.startswith(">"):
				seq = seq + Seq(line.strip("\n"))
	print len(seq)
	print seq.count("A")/float(len(seq))*100
	print seq.count("T")/float(len(seq))*100
	print seq.count("G")/float(len(seq))*100
	print seq.count("C")/float(len(seq))*100

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		get_ACGTcontent_from_FASTA(element)


