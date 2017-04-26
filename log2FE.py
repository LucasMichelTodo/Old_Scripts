#!/usr/bin/env python

# Import packages
import sys
import pandas
import numpy

def log2_transform(bedgraph):
	bdg = pandas.read_csv(filepath_or_buffer= bedgraph, sep="\t", header=None, index_col=None)
	bdg.loc[:,3] = bdg[3].apply(lambda x: numpy.log2(numpy.float64(x)))
	bdg.to_csv(path_or_buf= bedgraph+"_log2", sep="\t")
	print "Done!" 

if __name__ == "__main__":
	files = sys.argv[1:]
	for element in files:
		log2_transform(element)
