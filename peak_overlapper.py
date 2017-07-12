#!/usr/bin/env python

import argparse as ap
import pandas as pd


parser = ap.ArgumentParser(description='Reads peaks in PEAKS1.xls and keeps only those that overlap with a peak in PEAKS2.xls. Afterwards applies FE and QVAL filters.')
parser.add_argument('-p1',dest='peaks1', type=str, help='Peak file to de filtered.', required=True)
parser.add_argument('-p2',dest='peaks2', type=str, help='Peak file to be used as a filter.', required=True)
parser.add_argument('-f', '--fold-enrichment', dest='fe', type=float, help='Leave only peaks with fold-enrichment >= than FE (float, default: 1.5)', default=1.5)
parser.add_argument('-q', '--qval', dest='qval', type=float, help='Leave only peaks with -log10(qvalue) >= than QVAL (float, default: 15.0)', default=15.0)
parser.add_argument('--or-filters', dest='orfilter', help='Keep peaks passing either fe filter OR qval.', action='store_true', default=False)


if __name__ == "__main__":
	args = parser.parse_args()

# Load .xls files:
peaks1 = pd.read_csv(filepath_or_buffer= args.peaks1, sep="\t", header=0, index_col=None, comment="#")
peaks2 = pd.read_csv(filepath_or_buffer= args.peaks2, sep="\t", header=0, index_col=None, comment="#")

# Create a dictionary from contrast peaks file (file to be filtered). Each chromosome is a key with a list of ("start", "stop", "name") tupples representing peaks as value:
contrast = {}
for index, row in peaks1.iterrows():
	if row["chr"] in contrast.keys():
		contrast[row["chr"]].append((row["start"], row["end"], row["name"]))
	else: 
		contrast[row["chr"]] = [(row["start"], row["end"], row["name"])]

# Create a dictionary from reference peaks file (filter file):
reference = {}
for index, row in peaks2.iterrows():
	if row["chr"] in reference.keys():
		reference[row["chr"]].append((row["start"], row["end"]))
	else: 
		reference[row["chr"]] = [(row["start"], row["end"])]

# Add any keys to reference that may be missing as empty entries:
for key in contrast.keys():
	if key not in reference.keys():
		reference[key] = []

# Comparing the two dictionaries and keeping only overlapping peaks names (any partial overlap is kept):
hits = []
i = 0
for key in contrast:
	for element in contrast[key]:
		while i < len(reference[key]):
			if reference[key][i][0] <= element[1] and element[0] <= reference[key][i][1]:
				hits.append(element[2])
				i += 1
			else:
				i += 1
		i = 0

# Subsetting original peaks data.frame:
result = peaks1[peaks1["name"].isin(hits)]

# Appling FE and QVAL filters:
if args.orfilter:
	result = result[(result['fold_enrichment'] >= args.fe) | (result['-log10(qvalue)'] >= args.qval)]
else:
	result = result[(result['fold_enrichment'] >= args.fe) & (result['-log10(qvalue)'] >= args.qval)]

# Writing output to csv:
result.to_csv(path_or_buf=args.peaks1.replace(".xls", "_overlappandfilter.csv"), sep="\t", mode="w")



