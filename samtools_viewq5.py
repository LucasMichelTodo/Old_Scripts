#!/usr/bin/env python

import sys

if __name__ == "__main__":
	filenames = sys.argv[1:]
	print filenames
	for element in filenames:
		cmd = "samtools view -b -q 5 {} > {}" .format(element, element.replace(".bam", "_q5.bam"))
		subprocess.call(cmd, shell=True)
