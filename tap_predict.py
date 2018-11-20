#!/usr/bin/env python

import sys
import subprocess
from tqdm import tqdm

def tap_predict(fasta):

	cmd = ("python ~/Programs/netchop/predict.py --method netchop {}" .format(fasta) +
			" > {}_tap_predictions.txt" .format(fasta.replace(".fasta", "")))

	print cmd
	subprocess.call(cmd, shell=True)

if __name__ == "__main__":
	fileins = sys.argv[1:]
	for i in fileins:
		tap_predict(i)
