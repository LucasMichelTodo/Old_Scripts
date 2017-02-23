import sys
import subprocess
from tqdm import tqdm

filenames = sys.argv[1:]

for file in tqdm(filenames):
	cmd = "samtools view -h {} > {}" .format(file, file.replace(".bam", ".sam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools sort {} > {}" .format(file, file.replace(".bam", "_sorted.bam"))
	subprocess.call(cmd, shell=True)
	cmd = "samtools index {} > {}" .format(file, file.replace(".bam", "_sorted.bam.bai"))
	subprocess.call(cmd, shell=True)
