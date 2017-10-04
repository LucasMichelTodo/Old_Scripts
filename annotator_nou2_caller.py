#!/usr/bin/env python

# Import packages if __name__ == "__main__":
import sys
import subprocess as sp

filenames = sys.argv[1:]
print filenames
for element in filenames:
	cmd = "annotator_nou2.py {}" .format(element[2:])
	print cmd
	sp.call(cmd, shell=True)