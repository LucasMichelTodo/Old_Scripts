#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
import argparse as ap
import pysam as py
import pandas as pd
import numpy as np

parser = ap.ArgumentParser(description='Counting reads in peaks')
parser.add_argument('peaks', type=str, help='Peak files to be complemented.')

if __name__ == "__main__":
	args = parser.parse_args()

peaks1 = pd.read_csv(filepath_or_buffer= args.peaks, sep="\t", header=0, index_col=None, comment="#")
treat = "../../"+args.peaks.split("_")[0]+"_me_sort_q5.bam"
treat_bam = py.AlignmentFile(treat, "rb")
treat_max = []

for index,row in tqdm(peaks1.iterrows()):
	treat_cov = treat_bam.count_coverage(reference=row[0], start=row[1], end=row[2])
	treat_cov = np.sum(treat_cov, axis=0)
	treat_max.append(np.max(treat_cov))

treat_max.sort(reverse=True)

print treat_max

with open(args.peaks.split("_")[0]+"_maxpeaks.txt", "a+") as file:
	for element in treat_max:
		file.write("{}\n" .format(element))