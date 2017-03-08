#!/usr/bin/env python

import pandas

# Import data.frame with exons
exons = pandas.DataFrame.from_csv(path = "/home/lucas/ISGlobal/Gen_Referencies/true_exons.txt", sep=" ", header=None, index_col=None)
exons.columns=["chrom", "type", "start", "stop"]
grouped_df = exons.groupby("chrom")

# Group by chrmosome and sort
df = exons.sort_values(by="start",axis=0).groupby("chrom")

# Creating a dictionary with chromosomes as keys and start-stop positions as a list of tuples
coding = {}
for name, group in df:
	for idx, row in group.iterrows():
		if name in coding:
			coding[name].append((row["start"], row["stop"]))
		else:
			coding[name] = [(row["start"], row["stop"])]


# Purging ovelapping exons (alternative splicing). 
n_iter = [1,2] #setting number of needed iterations
for i in n_iter:
	for key in coding:
		for i in range(len(coding[key])-1):
				if coding[key][i][1] > coding[key][i+1][0]:
					coding[key][i] = (min([coding[key][i],coding[key][i+1]])[0], max([coding[key][i],coding[key][i+1]])[1])
					coding[key][i+1] = (0,0) # Set the ones to remove to (0,0)
				else:
					pass

for key in coding:
	coding[key] = [i for i in coding[key] if i[0] > 0] # Remove all (0,0)


# Checking no more overlaps are present.
for key in coding:
	for i in range(len(coding[key])-1):
			if coding[key][i][1] > coding[key][i+1][0]:
				print "Ojut!"
				print coding[key][i]
				print coding[key][i+1]
			else:
				pass

# Create non coding dictionary








