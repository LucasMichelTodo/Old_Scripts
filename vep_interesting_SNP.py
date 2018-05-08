#!/usr/bin/env python

# Import packages
import sys
import pybedtools as py
import subprocess as sp

def vep_snps(anotated_file):
    snps = []
    with open(anotated_file, "r+") as filein:
        for line in filein:
            # if line.startswith("\"Chrom\""):
            #     pass
            if line.startswith("Chrom"):
                pass
            else:
                # chrom = line.split(",")[0].replace("\"", "")
                # bp = line.split(",")[1]
                chrom = line.split()[0]
                bp = line.split()[1]
                snps.append((chrom,bp))

    ref_vcf = py.BedTool("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/All_SNP.vcf")

    for snp in snps:
        for variant in ref_vcf:
            if str(snp[0]) == str(variant.chrom) and str(snp[1]) == str(variant.start):
                print str(variant).strip()




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for element in filenames:
        vep_snps(element)
