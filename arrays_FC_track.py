#!/usr/bin/env python

import sys

gene_locations = {}
genes_FC = []

with open("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/whole_FC_df2.csv", "r+") as filein1:
    header = True
    for line in filein1:
        if header:
            header = False
        else:
            gene = line.strip().split(",")[3].replace("\"", "").split(": ")[1]
            try:
                fc1 = float(line.strip().split(",")[4])
            except:
                fc1 = 0.0
            try:
                fc2 = float(line.strip().split(",")[5])
            except:
                fc2 = 0.0
            try:
                fc3 = float(line.strip().split(",")[6])
            except:
                fc3 = 0.0

            FC = max([fc1, fc2, fc3], key=abs)

            variant = line.strip().split(",")[9]

            genes_FC.append((gene, FC, variant))


with open("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", "r+") as gff:
    for line in gff:
        if line.startswith("#"):
            pass
        elif line.split()[2] == "gene":

            gene = line.strip().split()[8].split(";")[0].replace("ID=", "")
            chrom = line.split()[0]
            start = int(line.split()[3])
            stop = int(line.split()[4])
            strand = line.split()[6]

            gene_locations[gene] = (chrom, start, stop, strand)


for entry in genes_FC:
    if entry[0] in gene_locations.keys():

        chrom = gene_locations[entry[0]][0]
        start = gene_locations[entry[0]][1]
        stop = gene_locations[entry[0]][2]
        strand = gene_locations[entry[0]][3]
        # fc1 = entry[1]
        # fc2 = entry[2]
        # fc3 = entry[3]
        FC = entry[1]
        variant = entry[2]

        # interval = (stop-start+1)/3

        # if (abs(fc1) >= 0.5849625) or (abs(fc2) >= 0.5849625) or (abs(fc3) >= 0.5849625):
        #     with open("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/Tracks/array_FC_track_alltimes.bed", "a+") as outfile:
        #         outfile.write("{}\t{}\t{}\t{}\n" . format(chrom, str(start), str(start+interval-1), fc1))
        #         outfile.write("{}\t{}\t{}\t{}\n" . format(chrom, str(start+interval), str(start+(2*interval)-1), fc2))
        #         outfile.write("{}\t{}\t{}\t{}\n" . format(chrom, str(start+(2*interval)), str(stop), fc3))

        if abs(FC) >= 0.5849625:
            with open("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/Tracks/array_FC_track_max.bed", "a+") as outfile:
                outfile.write("{}\t{}\t{}\t{}\n" . format(chrom, str(start), str(stop), FC))


        if variant == 'TRUE' and abs(FC) >= 0.5849625:
            with open("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/Tracks/array_FC_track_Variant.bed", "a+") as outfile:
                outfile.write("{}\t{}\t{}\t{}\n" . format(chrom, str(start), str(stop), FC))
