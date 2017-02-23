#!/usr/bin/env python

import sys
from tqdm import tqdm

filenames = sys.argv[1:]

def line_count(filename):
	with open(filename) as f: print sum(1 for line in tqdm(f))

for file in tqdm(filenames):
	line_count(file)