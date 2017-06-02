#!/usr/bin/env python

import sys
import pysam
import numpy

def get_factor(bam1, bam2):
	
	bams = {bam1:0, bam2:0}
	
	for key in bams:

		bam = pysam.AlignmentFile(key, "rb")
		coverage = []

		cov1 = bam.count_coverage(reference="Pf3D7_10_v3", start=1428000, end=1438000)
		cov2 = bam.count_coverage(reference="Pf3D7_12_v3", start=766000, end=786000)
		cov3 = bam.count_coverage(reference="Pf3D7_12_v3", start=903500, end=915600)
		cov4 = bam.count_coverage(reference="Pf3D7_06_v3", start=722500, end=745100)

		total_cov1 = numpy.sum(cov1, axis=0)
		total_cov2 = numpy.sum(cov2, axis=0)
		total_cov3 = numpy.sum(cov3, axis=0)
		total_cov4 = numpy.sum(cov4, axis=0)
 
		coverage = numpy.concatenate((total_cov1, total_cov2, total_cov3, total_cov4))
		bams[key] = numpy.mean(coverage)

	factor = bams[bam1]/bams[bam2]

	print bams
	print factor


if __name__ == "__main__":
	filenames = sys.argv[1:]
	get_factor(filenames[0], filenames[1])