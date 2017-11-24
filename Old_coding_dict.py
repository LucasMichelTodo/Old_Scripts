# Get reference annotated genome

ref_fasta = open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-30_Pfalciparum3D7_AnnotatedCDSs.fasta", "r+")
coding_regions = {}
for line in ref_fasta:
	if line.startswith(">"):
		gene = line.split("|")
		gene = map(lambda x: x.strip(" >\n"), gene)
		dict_string = map(lambda x: x.split("="), gene[1:])
		gene_dict = {}
		for element in dict_string:
			gene_dict[element[0]] = element[1]
		coding_regions[gene[0]] = gene_dict
coding_df = pandas.DataFrame.from_dict(data=coding_regions, orient = "index")


## Creating non-coding regions reference:

regions = coding_df["location"]

sizes = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
		else:
			pass

# sizes["Pf_M76611"] = sizes.pop("M76611")
# sizes["Pf3D7_API_v3"] = sizes.pop("PFC10_API_IRAB")

chrom = map(lambda x: x.split(":"), regions)
chroms = {}
for pair in chrom:
	if pair[0] in chroms:
		chroms[pair[0]].append(pair[1])
	else:
		chroms[pair[0]] = [pair[1]]

for key in chroms:
	chroms[key] = map(lambda x: x.strip("([+-])").split("-"), chroms[key])

noncoding = {}
for key in chroms:
	noncoding[key] = []
	for i in chroms[key]:
		noncoding[key].append(int(i[0])-1)
		noncoding[key].append(int(i[1])+1)

noncoding["M76611"] = noncoding.pop("Pf_M76611")
noncoding["PFC10_API_IRAB"] = noncoding.pop("Pf3D7_API_v3")

for key in noncoding:
	noncoding[key].insert(0,1)
	noncoding[key].append(sizes[key])

# Set one 0 to 1 (coding region starts at pos 1):
noncoding["PFC10_API_IRAB"][1] = 1
# noncoding["Pf3D7_API_v3"][1] = 1

## Remove overlapping genes:

n_iter = [0,1]
for i in n_iter:
	for key in noncoding:
		for i in range(len(noncoding[key])):
			if len(noncoding[key][i-2:i+2]) == 4:
				if noncoding[key][i] < noncoding[key][i-1]:
					# print "Overlap in Chrom: {} position {}" .format(key, i)
					# print noncoding[key][i-2:i+2]
					noncoding[key][i-2:i+2] = [min(noncoding[key][i-2:i+2]),max(noncoding[key][i-2:i+2])]
			else:
				pass

## Check no overlapping genes are still present:

for key in noncoding:
	for i in range(len(noncoding[key])):
		if noncoding[key][i] < noncoding[key][i-1] and i != 0 and i !=1:
			print "Overlap in Chrom: {} position {}" .format(key, i)
			print noncoding[key][i-2:i+2]