#!/usr/bin/env python

import argparse as ap
import pandas as pd


parser = ap.ArgumentParser(description='Adding information to Broad_Peaks.xls for further filtering')
parser.add_argument('-p1,',dest='peaks1', type=str, help='Peak files to de filtered.', required=True)
parser.add_argument('-p2',dest='peaks2', type=str, help='Peak files to be used as a filter.', required=True)

if __name__ == "__main__":
	args = parser.parse_args()

peaks1 = pd.read_csv(filepath_or_buffer= args.peaks1, sep="\t", header=0, index_col=None, comment="#")
peaks2 = pd.read_csv(filepath_or_buffer= args.peaks2, sep="\t", header=0, index_col=None, comment="#")

contrast = {}
for index, row in peaks1.iterrows():
	if row["chr"] in contrast.keys():
		contrast[row["chr"]].append((row["start"], row["end"], row["name"]))
	else: 
		contrast[row["chr"]] = [(row["start"], row["end"], row["name"])]

reference = {}
for index, row in peaks2.iterrows():
	if row["chr"] in reference.keys():
		reference[row["chr"]].append((row["start"], row["end"]))
	else: 
		reference[row["chr"]] = [(row["start"], row["end"])]

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


result = peaks1[peaks1["name"].isin(hits)]

result = result[(result['fold_enrichment'] >= 2) | (result['-log10(qvalue)'] >= 15)]

result.to_csv(path_or_buf=args.peaks1.replace(".xls", "_ovelappandfilter.csv"), sep="\t", mode="w")



