#!/usr/bin/env python

import pandas

exons = pandas.DataFrame.from_csv(path = "/home/lucas/ISGlobal/Gen_Referencies/true_exons.txt", sep=" ", header=None, index_col=None)
exons.columns=["chrom", "type", "start", "stop"]
grouped_df = exons.groupby("chrom")
for name, group in grouped_df:
	print name
	print group.sort_values(by="start",axis=0)



