#!/usr/bin/env python

from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo
from Bio import AlignIO
import sys
import pandas as pd
import scipy as sc
import scipy.stats

# Input a pandas series 
def ent(data): 
    p_data= data.value_counts()/len(data) # calculates the probabilities
    entropy=scipy.stats.entropy(p_data)  # input probabilities to get the entropy 
    return entropy


def most_common(lst):
    return max(set(lst), key=lst.count)


# Input a ClustalW .aln file
def ClustalEnt(clustal_file):
	aa_var = [] # container for list of aa in a position
	entropies = [] # container for final entropy result per position
	alignment = AlignIO.read(open(clustal_file), "clustal")
	
	with open(clustal_file.replace("_sorted.aln", "_consensus.fasta"), "w+") as file: # Creating header for consensus file. Setting "w+" ensures the file will be created from scratch if re-run.
		#for record in alignment:
		file.write(">"+alignment[0].id+"\n")

	#summary = SummaryInfo(alignment)
	#print summary.gap_consensus(threshold=0.4,ambiguous="*")[43]

	for i in range(0,alignment.get_alignment_length()):   
		for record in alignment:
			aa_var.append(record[i])
		data = pd.Series(aa_var)
		entropies.append(ent(data))
		
		#print aa_var
		#print most_common(aa_var)

		if ent(data) <= 0.5:
			 with open(clustal_file.replace("_sorted.aln", "_consensus.fasta"), "a+") as file:
		 		file.write(most_common(aa_var))
		else:
			with open(clustal_file.replace("_sorted.aln", "_consensus.fasta"), "a+") as file:
		 		file.write("*")

		aa_var = []

	#print entropies



if __name__ == "__main__":
	filenames= sys.argv[1:]
	print filenames
	for file in filenames:	
		ClustalEnt(file)
	
