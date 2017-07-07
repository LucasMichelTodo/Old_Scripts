#!/usr/bin/env python

import subprocess as sp

samples = ["C2", "A7K9", "E5K9", "10G", "1.2B"]

for sample in samples:
	cmd = "java -jar $PICARD SortSam INPUT={}_in_sort_q5.bam OUTPUT={}_GATKsorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" \
	.format(sample, sample)
	sp.call(cmd, shell=True)
	cmd = "java -jar $PICARD MarkDuplicates INPUT={}_GATKsorted.bam OUTPUT={}_GATKsorted_dupl.bam METRICS_FILE={}metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT" \
	.format(sample, sample, sample)
	sp.call(cmd, shell=True)
	cmd = "java -jar $PICARD AddOrReplaceReadGroups INPUT={}_GATKsorted_dupl.bam OUTPUT={}_GATKsorted_dupl_aorg.bam RGID=group{} RGLB=lib{} RGPL=illumina RGPU=unit1 RGSM={}_sample" \
	.format(sample, sample, sample, sample, sample)
	sp.call(cmd, shell=True)
	cmd = "java -jar $PICARD BuildBamIndex INPUT={}_GATKsorted_dupl_aorg.bam" \
	.format(sample)
	sp.call(cmd, shell=True)
	
# Merge tots aquests:
sp.call("java -jar $PICARD MergeSamFiles INPUT=C2_GATKsorted_dupl_aorg.bam INPUT=A7K9_GATKsorted_dupl_aorg.bam INPUT=E5K9_GATKsorted_dupl_aorg.bam INPUT=10G_GATKsorted_dupl_aorg.bam INPUT=1.2B_GATKsorted_dupl_aorg.bam OUTPUT=all_merged.bam SORT_ORDER=coordinate CREATE_INDEX=true", shell=True)
sp.call("java -jar $GATK -T RealignerTargetCreator -R ref.fasta -I all_merged.bam -o all_merged.bam.list", shell=True)
sp.call("java -jar $GATK -T IndelRealigner -R ref.fasta -I all_merged.bam -targetIntervals all_merged.bam.list -o all_merged_realigned.bam", shell=True)
sp.call("java -jar $PICARD FixMateInformation INPUT=all_merged_realigned.bam OUTPUT=all_merged_realigned_fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true", shell=True)
sp.call("java -jar $GATK -T  UnifiedGenotyper -R ref.fasta -I all_merged_realigned_fixed.bam -o All_SNP.vcf", shell=True)
