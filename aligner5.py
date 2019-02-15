#!/usr/bin/env python

import sys
import os
import subprocess

def align(filenames):

	#subprocess.call("mkdir Results_edited_genomes", shell=True)

	tupple_list = []
	i = 0
	while i < len(filenames):
		tupple_list.append((filenames[i], filenames[i+1]))
		i += 2

	for tupple in tupple_list:
		# cmd = "bbduk.sh in={} in2={} out1={} out2={} outm={} ref=/home/lucas/Programs/bbmap/resources/adapters.fa  ktrim=r k=22 mink=6 overwrite=t"\
		#  	.format(tupple[0], tupple[1], tupple[0].replace(".fastq", "_clean.fastq"), tupple[1].replace(".fastq", "_clean.fastq"), tupple[0].split("_")[0]+"_bad.fastq")
		# subprocess.call(cmd, shell=True)

		bw2_output = tupple[0].split("-")[0].split("_")[0]+"_multiplecopiestest.sam"

		cmd = "~/Programs/bowtie2-2.3.0-legacy/bowtie2 -p 4 \
		-x /media/lucas/Disc4T/Projects/Oriol/Plasmid_Multiple_Integrations/PlasmoDB-39_Pfalciparum3D7_Genome_3E_F1_Conca.fas\
		 -1 {} -2 {} > ../Plasmid/Aligns/{}"\
			.format(tupple[0],
					tupple[1],
					bw2_output)

		subprocess.call(cmd, shell=True)

		cmd = "samtools view -h -F 4 ../Plasmid/Aligns/{} > ../Plasmid/Aligns/mapped_{}" .format(bw2_output, bw2_output)
		subprocess.call(cmd, shell=True)

		cmd = "rm ../Plasmid/Aligns/{}" .format(bw2_output)
		subprocess.call(cmd, shell=True)

		cmd = "from_sam_to_bambai.py ../Plasmid/Aligns/mapped_{}" .format(bw2_output)
		subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    filenames = sys.argv[1:]
    align(filenames)
