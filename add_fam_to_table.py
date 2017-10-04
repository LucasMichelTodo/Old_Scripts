#!/usr/bin/env python

import csv
import sys
import pandas

fam = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Gen_Referencies/Families.csv", sep="\t", header=None, index_col=None)

def traduir(taula):
	fams = []
	with open(taula, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		rownum = 0
		for row in reader:
			if rownum != 0:
				fams.append(fam.loc[fam[2] == row[0]][1].values)
			rownum += 1


	print fams



if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		traduir(element)