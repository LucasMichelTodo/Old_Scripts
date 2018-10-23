#!/usr/bin/env python

import sys
import subprocess
from tqdm import tqdm

def predict_tcd4(fasta):

	## Edit this line to desired output directory
	output_dir = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/IEDB_recommended/Tcd8_predictions/"

	with open("/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/iedb_cow_alleles.txt", "r+") as infile:
		for line in tqdm(infile):
			hla = line.strip()

			cmd = ("python ~/Programs/mhc_i/src/predict_binding.py IEDB_recommended {} 9 {}" .format(hla, fasta) +
					" > {}{}_predictions.txt" .format(output_dir, hla.replace("/", "_")))

			print cmd
			subprocess.call(cmd, shell=True)

if __name__ == "__main__":
	filein = sys.argv[1]
	predict_tcd4(filein)
