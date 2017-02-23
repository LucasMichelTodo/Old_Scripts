#!/usr/bin/env python

import sys
import os
import subprocess
from tqdm import tqdm
from itertools import (takewhile,repeat)

filenames = sys.argv[1:]

for file in tqdm(filenames):
	if file[-4:] == ".sam":
		cmd = "head -220000 {} > {}" .format(file, file.replace(".sam", "_sample.sam"))
		subprocess.call(cmd, shell = True)
	elif file[-4:] == ".bam": 
		cmd = "head -220000 {} > {}" .format(file, file.replace(".bam", "_sample.bam"))
		subprocess.call(cmd, shell = True)
	else:
		print "Not a bam or sam file!"
