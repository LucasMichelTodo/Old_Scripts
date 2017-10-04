#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
import argparse as ap
import pysam as py
import pandas as pd
import numpy as np

	
depths = {"10G":4535787.0, "1.2B":4997472.0, "A7": 10915415.0, "C2": 5176848.0, "E5": 9366386.0}

parser = ap.ArgumentParser(description='Adding information to Broad_Peaks.xls for further filtering')
parser.add_argument('peaks', type=str, help='Peak files to be complemented.')
parser.add_argument('-t', type=str, help='Treatment bam file', required=True, dest="treat")
parser.add_argument('-c', type=str, help='Control bam file', required=True, dest="ctrl")

if __name__ == "__main__":
	args = parser.parse_args()

peaks1 = pd.read_csv(filepath_or_buffer= args.peaks, sep="\t", header=0, index_col=None, comment="#")
# treat = args.peaks.split("_")[0]+"_me_sort_q5.bam"
# ctrl = args.peaks.split("_")[1]+"_me_sort_q5.bam"

print "Treatment file: "+args.treat+"\n"
print "Control file: "+args.ctrl+"\n"

treat_bam = py.AlignmentFile(args.treat, "rb")
ctrl_bam = py.AlignmentFile(args.ctrl, "rb")

treat_max = []
treat_total = []
treat_mean = []

ctrl_max = []
ctrl_total = []
ctrl_mean = []

for index,row in tqdm(peaks1.iterrows()):
	#print "Chrom: {} \t Start: {} \t Stop: {}" .format(row[0], row[1], row[2])
	treat_cov = treat_bam.count_coverage(reference=row[0], start=row[1], end=row[2])
	treat_cov = np.sum(treat_cov, axis=0)
	ctrl_cov = ctrl_bam.count_coverage(reference=row[0], start=row[1], end=row[2])
	ctrl_cov = np.sum(ctrl_cov, axis=0)
	
	treat_max.append(np.max(treat_cov))
	treat_total.append(np.sum(treat_cov))
	treat_mean.append(np.mean(treat_cov))

	ctrl_max.append(np.max(ctrl_cov))
	ctrl_total.append(np.sum(ctrl_cov))
	ctrl_mean.append(np.mean(ctrl_cov))

for key in depths:
	if key in args.treat:
		treat_max = np.array(treat_max)/depths[key]*1000000
		treat_total = np.array(treat_total)/depths[key]*1000000
		treat_mean = np.array(treat_mean)/depths[key]*1000000
	if key in args.ctrl:
		ctrl_max = np.array(ctrl_max)/depths[key]*1000000
		ctrl_total = np.array(ctrl_total)/depths[key]*1000000
		ctrl_mean = np.array(ctrl_mean)/depths[key]*1000000

peaks1["Treat_max"] = treat_max
peaks1["Treat_total"] = treat_total
peaks1["Treat_mean"] = treat_mean
peaks1["Control_max"] = ctrl_max
peaks1["Control_total"] = ctrl_total
peaks1["Control_mean"] = ctrl_mean

peaks1.to_csv(path_or_buf=args.peaks.replace(".xls", "_filetered.csv"), sep="\t", mode="w")

