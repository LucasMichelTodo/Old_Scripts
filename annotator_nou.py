#!/usr/bin/env python

# Import packages 
import sys
import pybedtools as py

ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff")
ref = ref.sort()
gene_ref = ref.filter(lambda x: x[2] == "gene")



# prova = py.BedTool("/home/lucas/ISGlobal/prova.bed")

# result = prova.intersect(ref, wo=True)
# result2 = prova.closest(ref, D="ref")

# print prova.count()
# print result.filter(lambda x: x[12] == "gene").count()
# print result2.filter(lambda x: x[12] == "gene").count()

# print result2

def annotate_bed(bed_file):

	with open(bed_file.replace(".bed", "_annotated.csv"), "a+") as file:
		file.write("Gene\tGene-Coverage\t5-cov\t3-cov\tChrom\tAnnotations\n")

	bed = py.BedTool(bed_file)

	intersect = bed.intersect(gene_ref, wao=True)
	for line in intersect:
		if line[14] != str(-1):
			peak = {"gene":line[18], 
								"overlap":line[19],
								"gstart": int(line[13]),
								"gstop": int(line[14]),
								"start": int(line[1]),
								"stop" : int(line[2]),
								"gstrand" : line[16]}

			with open(bed_file.replace(".bed", "_annotated.csv"), "a+") as file:
				file.write(peak["gene"].split(";")[0].replace("ID=","")+"\t")
				file.write(peak["overlap"]+"\t")

				if peak["gstart"] - peak["start"] >= 0:
					if peak["gstrand"] == "+":
						tss = str(peak["gstart"] - peak["start"])
					else:
						tts = str(peak["gstart"] - peak["start"])

				if peak["stop"] - peak["gstop"] >= 0: 
					if peak["gstrand"] == "+":
						tts = str(peak["stop"] - peak["gstop"])
					else:
						tss = str(peak["stop"] - peak["gstop"])

				else:
					tss = str(0)
					tts = str(0)

				file.write(tss+"\t"+tts+"\t"+"\n")


		
		else:
			print "epa!"
			no_gene = py.BedTool(str(line), from_string=True)
			print no_gene.closest(ref, D="ref").filter(lambda x: x[22] == "gene")



if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		annotate_bed(element)
