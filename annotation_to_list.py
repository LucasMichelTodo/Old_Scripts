#!/usr/bin/env python

import csv
import sys

def from_anno_to_list(csv_file):

	idlist=[]
	with open(csv_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			for word in row:
				if "ID=" in word:
					idlist.append(word.replace("ID=", ""))


	with open(str(csv_file).replace("hits.csv", "list.txt"), "a+") as file:
		for element in idlist:
			file.write(element+"\n")



if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		from_anno_to_list(element)
