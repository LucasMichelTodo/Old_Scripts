#!/usr/bin/env python

from tqdm import tqdm
import re

in_pdb = re.compile("\d")
in_esmeraldo = re.compile("T")

clusters = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_05_03_18/clusters_70i.clstr", "r+") as infile:
	for line in infile:
		if line.startswith(">"):
			ID = line.strip()
			clusters[ID] = []
		else:
			clusters[ID].append(line.strip())

#print clusters

for key, value in tqdm(clusters.items()):
	proteomes = []
	for i in value:
		proteomes.append(i.split()[2][1:6])

	#print proteomes
	if len(set(proteomes)) < 7:
		del clusters[key]

print len(clusters)

# for key, value in tqdm(clusters.items()): #Iterate over clusters and return a string containing the first letter of each of their prot entries. The idea here is that all pdb entries start with a number while all esmeraldo entries start with a T.
# 	proteomes = ""
# 	for i in value:
# 		if i.split()[2].startswith(">"):
# 			proteomes += str(i.split()[2][1])
#
# 	if in_pdb.search(proteomes) and in_esmeraldo.search(proteomes):
# 		print proteomes, "MATCH!!"
# 	else:
# 		del clusters[key]
#
#
# print len(clusters)

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_05_03_18/clusters_70i_filtered.clstr", "a+") as outfile:
	for key, value in clusters.iteritems():
		outfile.write(key)
		outfile.write("\n")
		for i in value:
			outfile.write(i)
			outfile.write("\n")
