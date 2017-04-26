#!/usr/bin/env python

import sys
import os
import pandas
import re

gff_df = pandas.read_table(filepath_or_buffer = "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-30_Pfalciparum3D7.gff",\
 header = None, names=["Chrom", "Source", "Type", "Start", "Stop", ".", "Strand", "-", "Anno"], skiprows=18)

ID = re.compile("PF3D7_\d{7}")


if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		input_df = pandas.read_table(filepath_or_buffer = element)

		info = {}
		ls = []

		for line in input_df["Annotation"]:
			if ID.search(line):
				key = ID.search(line).group()
				query = re.compile(key)
				for row in gff_df["Anno"]:
					if query.search(row):
						description = re.compile(r"(?:\bdescription\b)(.*)")
						if description.search(row):
							info[key] = description.search(row).group(1)[1:]

		for line in input_df["Annotation"]:
			if ID.search(line):
				key = ID.search(line).group()
				ls.append(info[key])
			else:
				ls.append("No description.")
		
		#print info
		#print ls

		input_df["Description"] = ls

		df2 = input_df.rename(columns={input_df.columns.values[0]: 'Peak_ID'})
		final_df = df2.dropna(axis="columns", how="all")
		
		final_df.to_csv(path_or_buf=element.replace(".txt", "_curated.csv"), sep="\t")



