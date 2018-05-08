#!/usr/bin/env python

import sys
import pybedtools as py
from itertools import chain
from rosetta_to_dict import rosetta

gene_ref = py.BedTool("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff")

def anotate_variants(vcf_file):

    vcf = py.BedTool(vcf_file)
    intersect = vcf.intersect(gene_ref, wo=True) # Intersect each gene in the gff with the peaks in the bed file.

    #Print header before "for-loop":
    print "Chrom\tPos\tRef\tAlt\t1.2B_AD\t10G_AD\tA7_AD\tC2_AD\tE5_AD\t1.2B_RF\t10G_RF\tA7_RF\tC2_RF\tE5_RF\tGene\tRegion\tOld_fam\tAnnot"

    #Main loop:
    #Check if variant falls before, after or in ORF:
    for line in intersect:

        if int(line[1]) >= int(line[17])+int(line[23]) and int(line[1]) <= int(line[18])-int(line[24]):
            pos = "ORF"
        elif int(line[1]) < int(line[17])+int(line[23]):
            pos = "5region"
        elif int(line[1]) > int(line[18])-int(line[24]):
            pos = "3region"

        if line[19] == "-":
            if pos == "5region":
                pos = "3region"
            elif pos == "3region":
                pos = "5region"

        gene = line[22].split(";")[0].replace("ID=","")
        print_line = []
        print_line.append(line[0:2])
        print_line.append(line[3:5])

        #Get allelic depth from each sample. Must append inside a list so afterwards "chain" doesn't unlist content.
        #Calculate reference allele frequencies and put them on a list.
        rf_dp = []
        #Calculate depths and put them on a list (for filtering):
        for i in line[9:14]:
            if i == "./.":
                print_line.append([i])
                rf_dp.append((0,0))
                #rf.append(0)
                #dp.append(0)
            else:
                ad = i.split(":")[1]
                print_line.append([ad])
                try:
                    ad = map(int,ad.split(","))
                    rf_dp.append((round((float(ad[0])/sum(ad)),2),sum(ad)))
                    #rf.append(round((float(ad[0])/sum(ad)),2))
                    #dp.append(sum(ad))
                except:
                    rf_dp.append((0,0))
                    #rf.append(0)
                    #dp.append(0)

        for i in rf_dp:
            print_line.append([str(i[0])])

        print_line.append([gene,pos])

        try:
            print_line.append([rosetta[gene]["old_fam"]])
            print_line.append([rosetta[gene]["annot"]])
        except:
            print_line.append([".", "."])

        #Set fileters before printing the line:
        if any(x[0] < 0.3  and x [1] > 30 for x in rf_dp):

            #Calculate diferences in frequencie between each pair of diferent samples:
            difs = []
            for i in rf_dp:
                for j in rf_dp:
                    difs.append(abs(i[0]-j[0]))

            if any(x > 0.6 for x in set(difs)):

            # chain "unnlists" the list of lists we created.
                print "\t".join(list(chain.from_iterable(print_line)))
            else:
                pass



if __name__ == "__main__":
    filenames = sys.argv[1:]
    for element in filenames:
        anotate_variants(element)
