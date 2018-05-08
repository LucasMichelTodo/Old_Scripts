#!/usr/bin/env python

# Import packages
import sys
import pybedtools as py
import subprocess as sp
from itertools import chain
from rosetta_to_dict import rosetta

gene_ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff")

def annotate_from_genes(bed_file):

    bed = py.BedTool(bed_file) # Load bed file with peaks.
    intersect = bed.intersect(gene_ref, wo=True) # Intersect each gene in the gff with the peaks in the bed file.

    #Print header before "for-loop":
    print "Chrom\tPos\tRef\tAlt\t1_GT\t2_GT\t1_AD\t2_AD\t1_ratio\t2_ratio\t1_GQ\t2_GQ\tGene\tRegion\tOld_fam\tAnnot"
    #print "Chrom\tPos\tRef\tAlt\t1.2b_GT\t10g_GT\t1.2b_AD\t10g_AD\t1.2b_ratio\t10g_ratio\t1.2b_GQ\t10g_GQ\tGene\tRegion\tOld_fam\tAnnot"

    for line in intersect:

        if int(line[1]) >= int(line[16])+int(line[22]) and int(line[1]) <= int(line[17])-int(line[23]):
            pos = "ORF"
        elif int(line[1]) < int(line[16])+int(line[22]):
            pos = "5region"
        elif int(line[1]) > int(line[17])-int(line[23]):
            pos = "3region"

        if line[18] == "-":
            if pos == "5region":
                pos = "3region"
            elif pos == "3region":
                pos = "5region"

        gene = line[21].split(";")[0].replace("ID=","")
        print_line = []
        print_line.append(line[0:2])
        print_line.append(line[3:13])
        print_line.append([gene,pos])

        try:
            print_line.append([rosetta[gene]["old_fam"]])
            print_line.append([rosetta[gene]["annot"]])
        except:
            print_line.append([".", "."])

        # chain "unnlists" the list of lists we created.
        print "\t".join(list(chain.from_iterable(print_line)))



if __name__ == "__main__":
    filenames = sys.argv[1:]
    for element in filenames:
        annotate_from_genes(element)
