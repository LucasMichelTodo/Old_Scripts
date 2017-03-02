import sys
import os
from tqdm import tqdm
import subprocess

filenames = sys.argv[1:]
print filenames
for element in tqdm(filenames):
	name = element[:-4]
	cmd = "samtools view -bS {} > {}" .format(element, name+".bam")
	subprocess.call(cmd, shell = True)
	cmd = "samtools sort {} > {}" .format(name+".bam", name+"_sort.bam")
 	subprocess.call(cmd, shell = True)
 	cmd = "samtools index {} > {}" .format(name+"_sort.bam", name+"_sort.bam.bai")
 	subprocess.call(cmd, shell = True)
