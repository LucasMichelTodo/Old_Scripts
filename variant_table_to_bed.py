#!/usr/bin/env python

with open("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/a7_e5_Eli.csv", "r+") as file1:
    header = True
    for line in file1:
        if header:
            header = False
        else:
            line_list = line.strip().split("\t")
            line_list.insert(2, line_list[1])
            print "\t".join(line_list)
