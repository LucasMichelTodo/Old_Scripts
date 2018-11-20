#!/usr/bin/env python

import subprocess as sp

#samples = ["C2", "A7K9", "E5K9", "10G", "1.2B"]

#samples = ["3D7_me", "B11_me", "NF54_me"]

samples = ["3D7_ac", "B11_ac", "E5HA_ac", "NF54_ac"]

# inputPath = "/home/lucas/ISGlobal/Chip_Seq/Noves_dades/Results/Filtered_Bams_q5/"
#
# ### Pre-process each sample
# for sample in samples:
#
# 	# SortSam
# 	print "----------------Sorting SAMs----------------"
# 	cmd = "java -jar $PICARD SortSam INPUT={}{}_sort_q5.bam \
# 									OUTPUT={}_GATKsorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" .format(inputPath, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#MarkDuplicates
# 	print "----------------Marking duplicates----------------"
# 	cmd = "java -jar $PICARD MarkDuplicates INPUT={}_GATKsorted.bam OUTPUT={}_GATKsorted_dupl.bam METRICS_FILE={}metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT" \
# 	.format(sample, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#Erase pre-file:
# 	print "-----------------Erasing Pre-file: {}_GATKsorted.bam ----------------" .format(sample)
# 	cmd = "rm {}_GATKsorted.bam" .format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
#
# 	#AddOrReplaceReadGroups
# 	print "----------------Managing read-groups----------------"
# 	cmd = "java -jar $PICARD AddOrReplaceReadGroups INPUT={}_GATKsorted_dupl.bam OUTPUT={}_GATKsorted_dupl_aorg.bam RGID=group{} RGLB=lib{} RGPL=illumina RGPU=unit1 RGSM={}_sample" \
# 	.format(sample, sample, sample, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#Erase pre-file:
# 	print "-----------------Erasing Pre-file: {}_GATKsorted_dupl.bam ----------------" .format(sample)
# 	cmd = "rm {}_GATKsorted_dupl.bam" .format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#BuildBamIndex
# 	print "----------------Building index----------------"
# 	cmd = "java -jar $PICARD BuildBamIndex INPUT={}_GATKsorted_dupl_aorg.bam" \
# 	.format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#RealignerTargetCreator
# 	print "----------------Creating realignment list----------------"
# 	cmd = "java -jar $GATK -T RealignerTargetCreator -R ref.fasta -I {}_GATKsorted_dupl_aorg.bam -o {}_GATKsorted_dupl_aorg.bam.list" .format(sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#IndelRealigner
# 	print "----------------Realigning indels----------------"
# 	cmd = "java -jar $GATK -T IndelRealigner -R ref.fasta -I {}_GATKsorted_dupl_aorg.bam \
# 								-targetIntervals {}_GATKsorted_dupl_aorg.bam.list \
# 								-o {}_GATKsorted_dupl_aorg_realigned.bam" .format(sample, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#Erase pre-file:
# 	print "-----------------Erasing Pre-file: {}_GATKsorted_dupl_aorg.bam ----------------" .format(sample)
# 	cmd = "rm {}_GATKsorted_dupl_aorg.bam" .format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#FixMateInformation
# 	print "----------------Fixing mate information----------------"
# 	cmd = "java -jar $PICARD FixMateInformation INPUT={}_GATKsorted_dupl_aorg_realigned.bam \
# 												OUTPUT={}_GATKsorted_dupl_aorg_realigned_fixed.bam \
# 												SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" .format(sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#Erase pre-file:
# 	print "-----------------Erasing Pre-file: {}_GATKsorted_dupl_aorg_realigned.bam ----------------" .format(sample)
# 	cmd = "rm {}_GATKsorted_dupl_aorg_realigned.bam" .format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
#
# 	#### First Variant-calling round (run each sample separatedly)
# 	#UnifiedGenotyper
# 	print "\n----------------First Variant-Calling round----------------\n"
# 	cmd = "java -jar $GATK -T  UnifiedGenotyper -R ref.fasta -ploidy 2 -glm BOTH -I {}_GATKsorted_dupl_aorg_realigned_fixed.bam -o {}_SNP.vcf" .format(sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
#
# 	#### Create SNP and indels reference for second round of Variant_Calling
# 	# BaseRecalibrator
# 	print "\n----------------Creating VCF reference----------------\n"
# 	cmd = "java -jar $GATK -R ref.fasta -T BaseRecalibrator -I {}_GATKsorted_dupl_aorg_realigned_fixed.bam \
# 	 														-knownSites {}_SNP.vcf \
# 															-o {}_SNP_recal.grp" .format(sample, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
#
# 	#PrintReads
# 	print "\n----------------Creating new rearranged BAM----------------\n"
# 	cmd = "java -jar $GATK -R ref.fasta -T PrintReads -I {}_GATKsorted_dupl_aorg_realigned_fixed.bam \
# 														-BQSR {}_SNP_recal.grp \
# 														-o {}_SNP_recal_printreads_mq0.bam" .format(sample, sample, sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd
#
# 	#Erase pre-file:
# 	print "-----------------Erasing Pre-file: {}_GATKsorted_dupl_aorg_realigned_fixed.bam ----------------" .format(sample)
# 	cmd = "rm {}_GATKsorted_dupl_aorg_realigned_fixed.bam" .format(sample)
# 	sp.Popen(cmd, shell=True).wait()
# 	print cmd


# Modify this line so there are as many inputs as samples
### Merge all newly created bams:
# print "\n----------------Merging all rearranged samples----------------\n"
# cmd = "java -jar $PICARD MergeSamFiles INPUT={}_SNP_recal_printreads_mq0.bam \
# 									     INPUT={}_SNP_recal_printreads_mq0.bam \
# 										 INPUT={}_SNP_recal_printreads_mq0.bam \
# 										 INPUT={}_SNP_recal_printreads_mq0.bam \
# 										 OUTPUT=all_merged.bam SORT_ORDER=coordinate CREATE_INDEX=true" .format(*samples)
# sp.Popen(cmd, shell=True).wait()
# print cmd
#
# #Preprocess steps, again
# print "\n----------------Repeating pre-process steps with new merged BAM----------------\n"
# cmd = "java -jar $GATK -T RealignerTargetCreator -R ref.fasta -I all_merged.bam -o all_merged.bam.list"
# sp.Popen(cmd, shell=True).wait()
# print cmd
#
# cmd = "java -jar $GATK -T IndelRealigner -R ref.fasta -I all_merged.bam -targetIntervals all_merged.bam.list -o all_merged_realigned.bam"
# sp.Popen(cmd, shell=True).wait()
# print cmd
#
# #Erase pre-file:
# print "-----------------Erasing Pre-file: all_merged.bam ----------------"
# cmd = "rm all_merged.bam"
# sp.Popen(cmd, shell=True).wait()
# print cmd



cmd = "java -jar $PICARD FixMateInformation INPUT=all_merged_realigned.bam OUTPUT=all_merged_realigned_fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
sp.Popen(cmd, shell=True).wait()
print cmd

# #Erase pre-file:
# print "-----------------Erasing Pre-file: all_merged_realigned.bam ----------------"
# cmd = "rm all_merged_realigned.bam"
# sp.Popen(cmd, shell=True).wait()
# print cmd

#### Final Variant-Calling
print "\n----------------Final Variant-Calling----------------\n"
cmd = "java -jar $GATK -T  UnifiedGenotyper -R ref.fasta -ploidy 2 -glm BOTH \
												-stand_call_conf 0 -mbq 0 --min_quality_score 0 --max_quality_score 100 \
												-I all_merged_realigned_fixed.bam -o All_SNP.vcf"
sp.Popen(cmd, shell=True).wait()
print cmd
