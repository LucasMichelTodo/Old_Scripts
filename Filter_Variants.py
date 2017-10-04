#!/usr/bin/env python

import pandas as pd

table = pd.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/filtered_results.table", sep="\t")

table["1.2B_sample.AD"] = table["1.2B_sample.AD"].astype(str)
table["10G_sample.AD"] = table["10G_sample.AD"].astype(str)
table["A7K9_sample.AD"] = table["A7K9_sample.AD"].astype(str)
table["E5K9_sample.AD"] = table["E5K9_sample.AD"].astype(str)
table["C2_sample.AD"] = table["C2_sample.AD"].astype(str)

mask12B = []
for index, row in table.iterrows():
	if row["1.2B_sample.AD"] != "nan":
		if int(row["1.2B_sample.AD"].split(",")[0]) < int(row["1.2B_sample.AD"].split(",")[1]) and int(row["1.2B_sample.AD"].split(",")[0]) + int(row["1.2B_sample.AD"].split(",")[1]) > 20:
			mask12B.append(True)
		else:
			mask12B.append(False)
	else:
		mask12B.append(False)

mask10G = []
for index, row in table.iterrows():
	if row["10G_sample.AD"] != "nan":
		if int(row["10G_sample.AD"].split(",")[0]) < int(row["10G_sample.AD"].split(",")[1]) and int(row["10G_sample.AD"].split(",")[0]) + int(row["10G_sample.AD"].split(",")[1]) > 20:
			mask10G.append(True)
		else:
			mask10G.append(False)
	else:
		mask10G.append(False)

maskA7 = []
for index, row in table.iterrows():
	if row["A7K9_sample.AD"] != "nan":
		if int(row["A7K9_sample.AD"].split(",")[0]) < int(row["A7K9_sample.AD"].split(",")[1]) and int(row["A7K9_sample.AD"].split(",")[0]) + int(row["A7K9_sample.AD"].split(",")[1]) > 20:
			maskA7.append(True)
		else:
			maskA7.append(False)
	else:
		maskA7.append(False)

maskC2 = []
for index, row in table.iterrows():
	if row["C2_sample.AD"] != "nan":
		if int(row["C2_sample.AD"].split(",")[0]) < int(row["C2_sample.AD"].split(",")[1]) and int(row["C2_sample.AD"].split(",")[0]) + int(row["C2_sample.AD"].split(",")[1]) > 20:
			maskC2.append(True)
		else:
			maskC2.append(False)
	else:
		maskC2.append(False)

maskE5 = []
for index, row in table.iterrows():
	if row["E5K9_sample.AD"] != "nan":
		if int(row["E5K9_sample.AD"].split(",")[0]) < int(row["E5K9_sample.AD"].split(",")[1]) and int(row["E5K9_sample.AD"].split(",")[0]) + int(row["E5K9_sample.AD"].split(",")[1]) > 20:
			maskE5.append(True)
		else:
			maskE5.append(False)
	else:
		maskE5.append(False)


result = table[mask12B or mask10G or maskA7 or maskC2 or maskE5]
result.to_csv(path_or_buf="/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/custom_filtered_results.table", sep="\t", mode="w+")
