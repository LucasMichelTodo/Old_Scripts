#!/usr/bin/env python

import sys

reference = {}

with open("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/very_old_array_description.csv", "a+") as old_array:
    with open("/home/lucas/ISGlobal/Arrays/new_array_description_final.csv", "r+") as new_array:
        for line in new_array:
            linelist = line.strip().split("\t")
            reference[linelist[2].replace("\"","")] = (linelist[8],linelist[9])

        for line in old_array:
            linelist = line.strip().split("\t")

            if linelist[1] in reference.keys():
                print line.strip()+"\t"+reference[linelist[1]][0]+"\t"+reference[linelist[1]][1]

            else:
                print line.strip()+"\t"+linelist[2]+"\tNA"
