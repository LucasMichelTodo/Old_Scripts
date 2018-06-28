#!/usr/bin/env python

from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo
from Bio import AlignIO
import sys
import pandas as pd
import scipy as sc
import scipy.stats

# Calculate entropy given a panda's series object
def ent(data):
    p_data= data.value_counts()/len(data) # calculates the probabilities
    entropy=scipy.stats.entropy(p_data)  # input probabilities to get the entropy
    return entropy

# Get most common element from a list
def most_common(lst):
    return max(set(lst), key=lst.count)

# Calculate lenvenstein distance between two strings
def lv_number(a, b):
    string1 = a
    string2 = b
    distance = 0
    n1 = len(string1)
    n2 = len(string2)

    if n1 >= n2:
        for i in range(n1):
            if i < n2:
                if string1[i] != string2[i]:
                    distance += 1
            else:
                distance += 1
    else:
        for i in range(n2):
            if i < n1:
                if string2[i] != string1[i]:
                    distance -= 1
            else:
                distance -= 1


    return distance


# Input a ClustalW .aln file
def ClustalEnt(clustal_files):

    epitopes = []

    for clustal in clustal_files:

        aa_var = [] # container for list of aa in a position
        entropies = [] # container for final entropy result per position
        consensus = [] # container for consenses sequence
        alignment = AlignIO.read(open(clustal), "clustal")

        # For each position in the alignment create a list with all aa present in that position and calculate entropy over it.
        for i in range(0,alignment.get_alignment_length()):
            for record in alignment:
                aa_var.append(record.seq[i])
            data = pd.Series(aa_var)
            entropies.append(ent(data))

            # Create list with consensus aa for each position
            consensus.append(most_common(aa_var))

            # Restart container
            aa_var = []

        consensus_str = "".join(consensus)

        # Calculate which of the sequences in the alignment is more similar tot he consensus one.
        current_min = 99999999999999999999999999999999999
        current_seq = ""
        current_ref = ""

        # Define scores to prioritize some sources in case of equal distance to consensus.
        ref_order = {'TcCLB':5, 'AODP0':4, 'AYLP0':7, 'AQHO0':2, 'TcChr':8, 'Tcruz':6, 'ANOX0':1, 'TcX10':3, 'Expos':0, 'TcMAR':9}

        # Loop over all sequences in the alignment and minimize lenvenstein distance to consensus.
        for record in alignment:
            if lv_number(str(record.seq), consensus_str) < current_min:
                current_min = lv_number(str(record.seq), consensus_str)
                current_seq = record.seq
                current_ref = record.id
            elif lv_number(str(record.seq), consensus_str) == current_min:
                if ref_order[record.id[0:5]] < ref_order[current_ref[0:5]]:
                    current_min = lv_number(str(record.seq), consensus_str)
                    current_seq = record.seq
                    current_ref = record.id

        # Merge sequence and entropy information in a sigle (sliceable) array.
        final_seq_ent = []
        for i in range(0, len(current_seq)):
            final_seq_ent.append([current_seq[i], entropies[i]])

        # Break sequences by gaps, and only consider resulting sequences if len > 15.
        fragments = []
        fragment = []

        for i in final_seq_ent:
            if i[0] not in "-*":
                fragment += [i]
            else:
                if len(fragment) >= 15:
                    fragments.append(fragment)
                fragment = []
        if len(fragment) >= 15:
            fragments.append(fragment)

        #Split resulting fragments into 15bp epitopes and calculate global entropy (as the sum of each aa entropy in the epitope).
        for fragment in fragments:
            start = 0
            while start+15 <= len(fragment):
                epitope = "".join([i[0] for i in fragment[start:start+15]])
                entropy = sum([i[1] for i in fragment[start:start+15]])
                epitopes.append([epitope, entropy, current_ref])
                start += 1

    # Sort epitopes by entropy:
    epitopes = sorted(epitopes, key=lambda x: x[1])

    # Create references for B epitopes:
    epitopes_b = []
    with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/materias_primas/bcell_epitopes_to_cross.txt", "r+") as file:
    	for line in file:
    		epi = str(line.strip())
    		if len(epi) > 1:
    			epitopes_b.append(epi)


    for i in epitopes:
        if i[0] in epitopes_b:
            print i



if __name__ == "__main__":
    filenames= sys.argv[1:]
    ClustalEnt(filenames)
