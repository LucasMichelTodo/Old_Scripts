#!/usr/bin/env python

import sys

vep_anot = {}

def merge_vep_annotation(variants_file, vep_anotated_file):
    with open(vep_anotated_file, "r+") as filein1:
        for line in filein1:
            if line.startswith("#"):
                pass
            elif line.split()[10] != "-":
                chrom = line.split()[1].split(":")[0]
                pos = line.split()[1].split(":")[1].split("-")[0]
                var_type = line.split()[6]
                aa = line.split()[10]
                nuc = line.split()[11]
                vep_anot[(chrom, pos)] = (var_type, aa, nuc)

    with open(variants_file, "r+") as filein2:
        for line in filein2:
            chrom = line.split()[0]
            pos = line.split()[1]

            if line.startswith("Chrom"):
                print line.strip()+"\tVar_Type\tAAs\tCodons"
            else:
                if (chrom, pos) in vep_anot.keys():
                    anot = vep_anot[(chrom, pos)]
                    print line.strip()+"\t{}\t{}\t{}" .format(anot[0], anot[1], anot[2])
                else:
                    print line.strip()+"\t.\t.\t."




if __name__ == "__main__":
    filenames = sys.argv[1:]
    merge_vep_annotation(filenames[0], filenames[1])
