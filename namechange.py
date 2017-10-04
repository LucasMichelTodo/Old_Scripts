#!/usr/bin/env python

import sys
import os
import subprocess

filenames = sys.argv[1:]

lis = []
for file in filenames:
	 lis.append(file.split("-")[1])

lis2 = []
for file in lis:
	lis2.append(file.split("_"))

lis3 = []
for element in lis2:
	lis3.append(element[0][0:-2]+"_"+element[0][-2:]+"_"+element[3])

for i in range(len(filenames)):
	cmd = "mv {} ./{}" .format(filenames[i], lis3[i])
	subprocess.call(cmd, shell=True)

