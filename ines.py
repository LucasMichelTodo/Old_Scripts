#!/usr/bin/env python

#with open("/home/lucas/ISGlobal/tangolinux_split.txt", "r+") as infile1:
with open("/home/lucas/Programs/ines_2.bat", "r+") as infile1:
    for line in infile1:
        line_list = line.split()
        # name = line_list[1]
        # rename = name.split("|")[1]
        # line_list[1] = rename
        # if line_list[2] != "nt=\"N\"":
        #     del line_list[2]

        if line_list[2] != "nt=\"N\"":
            print line
        #print " ".join(line_list)
