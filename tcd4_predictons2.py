#!/usr/bin/env python

import sys
import subprocess
from tqdm import tqdm

def predict_tcd4(fasta):

	## Edit this line to desired output directory
	output_dir = "/media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Non_exposed/Tcd4_predictions/"

	with open("/media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Non_exposed/hla_ref_set.class_ii.txt", "r+") as infile:
		for line in tqdm(infile):
			hla = line.split("\t")[0].strip()

			cmd = ("python ~/Programs/mhc_ii/mhc_II_binding.py IEDB_recommended {} {}" .format(hla, fasta) +
					" > {}{}_predictions.txt" .format(output_dir, hla.replace("/", "_")))

			print cmd
			subprocess.call(cmd, shell=True)

if __name__ == "__main__":
	filein = sys.argv[1]
	predict_tcd4(filein)
