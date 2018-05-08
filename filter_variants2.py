#!/usr/bin/env python

import pandas as pd

ref = []
header = True
pd.set_option('display.expand_frame_repr', False)

## Previament hem de parsejar l'arxiu vcf!
## java -jar $GATK -R ref.fasta -T VariantsToTable -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF GQ -V haplotypecaller_SNP_ploidy.vcf -o haplotypecalle_SNP_ploidy_table.txt

## Obrir arxiu vcf (resultat del variant calling) i parsejar totes les entrades coma diccionaris en una llista.
with open("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/all_SNP_eli_table.txt", "r+") as file1:
		for line in file1:
			if header:
				header = False
			else:
				ref.append({"a7_e5":{"Chrom":line.split()[0],
								"Pos":line.split()[1],
								"Ref":line.split()[2],
								"Alt":line.split()[3],
								"A7_GT":line.split()[10],
								"E5_GT":line.split()[16],
								"A7_AD":line.split()[11],
								"E5_AD":line.split()[17],
								"A7_GQ":line.split()[12],
								"E5_GQ":line.split()[18]},
							"g10_b12":{"Chrom":line.split()[0],
									"Pos":line.split()[1],
									"Ref":line.split()[2],
									"Alt":line.split()[3],
									"1.2B_GT":line.split()[4],
									"10G_GT":line.split()[7],
									"1.2B_AD":line.split()[5],
									"10G_AD":line.split()[8],
									"1.2B_GQ":line.split()[6],
									"10G_GQ":line.split()[9]}})

filter_na = [item for item in ref if item["a7_e5"]["A7_AD"] != "NA" and item["a7_e5"]["E5_AD"] != "NA"]
filter_na = [item for item in filter_na if item["g10_b12"]["1.2B_AD"] != "NA" and item["g10_b12"]["10G_AD"] != "NA"]

a7_e5 = []
for i in filter_na:
	if (int(i["a7_e5"]["A7_AD"].split(",")[1]) + int(i["a7_e5"]["A7_AD"].split(",")[0]) > 10 and
		int(i["a7_e5"]["E5_AD"].split(",")[1]) + int(i["a7_e5"]["E5_AD"].split(",")[0]) > 10 and
		int(i["a7_e5"]["A7_GQ"]) >= 20 and int(i["a7_e5"]["E5_GQ"]) >= 20):
		a7_e5.append(i["a7_e5"])

for i in a7_e5:
	i["A7_ratio"] = float(i["A7_AD"].split(",")[1])/(int(i["A7_AD"].split(",")[0])+ float(i["A7_AD"].split(",")[1]))
	i["E5_ratio"] = float(i["E5_AD"].split(",")[1])/(int(i["E5_AD"].split(",")[0]) + float(i["E5_AD"].split(",")[1]))

dif_GT = [i for i in a7_e5 if i["A7_GT"] != i["E5_GT"]]
dif_a7_e5 = [i for i in a7_e5 if abs(i["A7_ratio"]-i["E5_ratio"]) > 0.3 or i["A7_ratio"] > 0.3 or i["E5_ratio"] > 0.3] #Edit this line to change diference bwtween strain ratios to filter.

df = pd.DataFrame(dif_a7_e5)
#df = pd.DataFrame(dif_GT)
#df = pd.DataFrame(a7_e5)
a7_e5_ordered = df[["Chrom", "Pos", "Ref", "Alt", "A7_GT", "E5_GT", "A7_AD", "E5_AD", "A7_ratio", "E5_ratio", "A7_GQ", "E5_GQ"]]
print a7_e5_ordered
a7_e5_ordered.to_csv(path_or_buf="/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/a7_e5.csv", sep="\t", mode="w+", index = False)


g10_b12 = []
for i in filter_na:
	if (int(i["g10_b12"]["1.2B_AD"].split(",")[1]) + int(i["g10_b12"]["1.2B_AD"].split(",")[0]) > 10 and
		int(i["g10_b12"]["10G_AD"].split(",")[1]) + int(i["g10_b12"]["10G_AD"].split(",")[0]) > 10 and
		int(i["g10_b12"]["1.2B_GQ"]) >= 20 and int(i["g10_b12"]["10G_GQ"]) >= 20):
		g10_b12.append(i["g10_b12"])

for i in g10_b12:
	i["1.2B_ratio"] = float(i["1.2B_AD"].split(",")[1])/(int(i["1.2B_AD"].split(",")[0])+ float(i["1.2B_AD"].split(",")[1]))
	i["10G_ratio"] = float(i["10G_AD"].split(",")[1])/(int(i["10G_AD"].split(",")[0]) + float(i["10G_AD"].split(",")[1]))

dif_12b_10g_GT = [i for i in g10_b12 if i["1.2B_GT"] != i["10G_GT"]]
dif_12b_10g = [i for i in g10_b12 if abs(i["1.2B_ratio"]-i["10G_ratio"]) > 0.3 or i["1.2B_ratio"] > 0.3 or i["10G_ratio"] > 0.3] #Edit this line to change diference bwtween strain ratios to filter.

df_12b_10g = pd.DataFrame(dif_12b_10g)
#df_12b_10g = pd.DataFrame(dif_12b_10g_GT)
#df_12b_10g = pd.DataFrame(g10_b12)
b12_g10_ordered = df_12b_10g[["Chrom", "Pos", "Ref", "Alt", "1.2B_GT", "10G_GT", "1.2B_AD", "10G_AD", "1.2B_ratio", "10G_ratio", "1.2B_GQ", "10G_GQ"]]
print b12_g10_ordered

b12_g10_ordered.to_csv(path_or_buf="/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/1.2b_10g.csv", sep="\t", mode="w+", index = False)
