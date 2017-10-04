#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import packages 
import sys
import pybedtools as py
import subprocess as sp
import pandas as pd

ref = pd.read_csv("/home/lucas/ISGlobal/prova_per_allargar_fins_seguent_gen.gff", sep='\t', header=None)
ref.columns = ["gene", "source", "feature", "start", "stop", "score", "strand", "frame", "attributes"]

print ref


# def valuation_formula(x, y):
#     return x * y * 0.5

# df['price'] = df.apply(lambda row: valuation_formula(row['x'], row['y']), axis=1)

# for i, row in df.iterrows():
#   ifor_val = something
#   if <condition>:
#     ifor_val = something_else
#   df.set_value(i,'ifor',ifor_val)

for i, row in ref.iterrows():
	if i < 23:
		elongation = 1000
		if ref.iloc[i+1]["stop"] - ref.iloc[i]["start"] < 1000:
			elongation = ref.iloc[i+1]["start"] - ref.iloc[i]["stop"]
		ref.set_value(i, "stop", ref.iloc[i]["stop"]+1000)

print ref



# ref2 = open("/home/lucas/ISGlobal/prova_per_allargar_fins_seguent_gen.gff", "r")

# current_gene = "Pf3D7_01_v3"

# for line in ref2:
# 	gene = line.split("\t")[0]
# 	start = line.split("\t")[3]
# 	stop = line.split("\t")[4]



