#!/usr/bin/env python

import re

fams = {}

with open("/home/lucas/ISGlobal/Gen_Referencies/variant_families.txt", "r+") as file1:	
	for line in file1:
		regex = re.compile(str(line.strip()))
		fams[(regex, regex.pattern)] = ""

with open("/home/lucas/ISGlobal/Gen_Referencies/families_all.txt", "r+") as file2:
	for line in file2:
		for key in fams:
			m = key[0].match(line)
			if m:
				print key[1]+"--->"+str(line)






