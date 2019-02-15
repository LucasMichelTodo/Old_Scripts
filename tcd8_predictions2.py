#!/usr/bin/env python

import sys
import subprocess
from tqdm import tqdm

def predict_tcd8(fasta):

	## Edit this line to desired output directory
	output_dir = "/media/lucas/Disc4T/Projects/tcruzi_Actual/Run_Non_exposed/Tcd8_predictions/"

	with open("/media/lucas/Disc4T/Projects/tcruzi_Actual/hla_ref.csv", "r+") as infile:
		for line in tqdm(infile):
			hla = line.split("\t")[0].strip()

			cmd = ("python ~/Programs/mhc_i/src/predict_binding.py IEDB_recommended {} 9 {}" .format(hla, fasta) +
					" > {}{}_predictions.txt" .format(output_dir, hla.replace("/", "_")))

			print cmd
			subprocess.call(cmd, shell=True)

if __name__ == "__main__":
	filein = sys.argv[1]
	predict_tcd8(filein)
