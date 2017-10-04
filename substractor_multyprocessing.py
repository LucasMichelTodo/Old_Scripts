#!/usr/bin/env python

import pysam
import csv 
import sys
import numpy as np
import pandas as pd

import threading
import os
import time
import multiprocessing

def coverage(bamfile):

	sam = pysam.AlignmentFile(bamfile, "rb")

	chrom_cov = {}
	sizes = {}
	
	with open("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes", "rb") as csvfile:
		chrom_sizes = csv.reader(csvfile, delimiter='\t')
		for row in chrom_sizes:
			if len(row) == 2:
				sizes[row[0]] = row[1]
			else:
				pass

	for key in sizes:
		coverage = []
		cov = sam.count_coverage(reference=key, start=1, end=int(sizes[key]))
		total_cov = numpy.sum(cov, axis=0)
		print "total_cov"
		print total_cov
		smooth_cov = numpy.array(pd.Series(total_cov).rolling(2000, min_periods=1).mean().tolist())
		chrom_cov[key] = smooth_cov
		print "smooth_cov"
		print smooth_cov

	return chrom_cov



NUM_WORKERS = 4

processes = [multiprocessing.Process(target=crunch_numbers) for _ in range(NUM_WORKERS)]
[process.start() for process in processes]
[process.join() for process in processes]
