#!/usr/bin/env python

import pandas as pd

ref = []
header = True
pd.set_option('display.expand_frame_repr', False)

## Previament hem de parsejar l'arxiu vcf!
## java -jar $GATK -R ref.fasta -T VariantsToTable -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF GQ -V haplotypecaller_SNP_ploidy.vcf -o haplotypecalle_SNP_ploidy_table.txt

## Obrir arxiu vcf (resultat del variant calling) i parsejar totes les entrades coma diccionaris en una llista.
with open("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/all_recal_table.txt", "r+") as file1:
		for line in file1:
			if header:
				header = False
			else:
				ref.append({"Chrom":line.split()[0],
							"Pos":line.split()[1],
							"Ref":line.split()[2],
							"Alt":line.split()[3],
							"A7_GT":line.split()[10],
							"E5_GT":line.split()[16],
							"A7_AD":line.split()[11],
							"E5_AD":line.split()[17],
							"A7_GQ":line.split()[12],
							"E5_GQ":line.split()[18],})

filter_na = [item for item in ref if item["A7_AD"] != "NA" and item["E5_AD"] != "NA"]

a7_e5 = []
for i in filter_na:
	if ((int(i["A7_AD"].split(",")[1])/(int(i["A7_AD"].split(",")[0])+1.0) >= 0 or
		int(i["E5_AD"].split(",")[1])/(int(i["E5_AD"].split(",")[0])+1.0) >= 0) and
		int(i["A7_AD"].split(",")[1]) + int(i["A7_AD"].split(",")[0]) > 10 and
		int(i["E5_AD"].split(",")[1]) + int(i["E5_AD"].split(",")[0]) > 10 and
		int(i["A7_GQ"]) >= 20 and int(i["E5_GQ"]) >= 20):
		a7_e5.append(i)

for i in a7_e5:
	i["A7_ratio"] = float(i["A7_AD"].split(",")[1])/(int(i["A7_AD"].split(",")[0])+ float(i["A7_AD"].split(",")[1]))
	i["E5_ratio"] = float(i["E5_AD"].split(",")[1])/(int(i["E5_AD"].split(",")[0]) + float(i["E5_AD"].split(",")[1]))

dif_GT = [i for i in a7_e5 if i["A7_GT"] != i["E5_GT"]]
dif_a7_e5 = [i for i in a7_e5 if abs(i["A7_ratio"]-i["E5_ratio"]) > 0.3] #Edit this line to change diference bwtween strain ratios to filter.

df = pd.DataFrame(dif_a7_e5)
#df = pd.DataFrame(dif_GT)
#df = pd.DataFrame(a7_e5)
df_ordered = df[["Chrom", "Pos", "Ref", "A7_GT", "E5_GT", "A7_AD", "E5_AD", "A7_ratio", "E5_ratio", "A7_GQ", "E5_GQ"]]
print df_ordered


df_ordered.to_csv(path_or_buf="/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/a7_e5_all_recal_table.txt", sep="\t", mode="w+", index = False)
