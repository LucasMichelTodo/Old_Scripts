#!/usr/bin/env python

from tqdm import tqdm

clusters = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_18_01_18/tcruzi_filetered_unique.clstr", "r+") as infile:
	for line in infile:
		if line.startswith(">"):
			ID = line.strip()
			clusters[ID] = []
		else:
			clusters[ID].append(line.strip())




# for key, value in tqdm(clusters.items()):
# 	proteomes = []
# 	for i in value:
# 		proteomes.append(i.split()[2][1:3])
# 	if len(set(proteomes)) < 2:
# 		del clusters[key]


print len(clusters)

# with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_proteomes/tcruzi_proteome_unique.clstr", "a+") as outfile:
# 	for key, value in clusters.iteritems():
# 		outfile.write(key)
# 		outfile.write("\n")
# 		for i in value:
# 			outfile.write(i)
# 			outfile.write("\n")

