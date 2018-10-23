#!/usr/bin/env python

import sys

def read_csv_and_cat(csv_file):
    col1 = []
    col2 = []
    with open(casv_file) as infile:
        for line in infile:
            col1.append(line.strip().split(",")[0])
            col2.append(line.strip().split(",")[1])

        with open("all_transposed_csv.csv", "a+") as outfile:
            outfile.write(",".join(col1)+"\n")
            outfile.write(",".join(col2)+"\n")

if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		read_csv_and_cat(file)







¡
