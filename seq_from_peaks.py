#!/usr/bin/env python

import sys
import csv

fasta_ref = {}
chrom = ""
seq = ""
with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Genome.fasta", "r+") as file1:	
	for line in file1:
		if line.startswith(">"):
			if chrom:
				fasta_ref[chrom] = seq
				seq = ""
			chrom = line.split(" | ")[0].replace(">", "")			
		else:
			if seq:
				seq = seq+line.strip()
			else:
				seq = line.strip()


for key, value in fasta_ref.iteritems():
	print key
	print value[0:10]

sizes = {}
with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
	chrom_sizes = csv.reader(csvfile, delimiter='\t')
	for row in chrom_sizes:
		if len(row) == 2:
			sizes[row[0]] = row[1]
		else:
			pass

print sizes

def get_2000_prior_bp(peaks_file):
	peaks = []
	with open(peaks_file, "r+") as file:
		for line in file:
			if line.startswith("#") or line.startswith("chr"):
				pass
			elif line.split():
				chrom = line.split()[0]

				if int(line.split()[1])-2000 > 0:
					start = int(line.split()[1])-2000
				else: 
					start = 1

				if int(line.split()[1])+2000 < sizes[chrom]:
					stop = int(line.split()[1])+2000
				else: 
					stop = sizes[chrom]

				peaks.append((chrom, start, stop))



	ext = "." + str(peaks_file.split(".")[1])
	i = 1
	for element in peaks:
		with open(peaks_file.replace(ext, "_prepeaks_seqs.fasta"), "a+") as outfile:
			seq = ">peak_"+str(i)
			outfile.write(seq)
			outfile.write("\n")
			outfile.write(fasta_ref[element[0]][element[1]-1:element[2]-1])
			outfile.write("\n")
			i += 1




if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		get_2000_prior_bp(element)