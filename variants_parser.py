#!/usr/bin/env python

import sys
from itertools import compress

def parse_variants(vars_table):

    header = True
    with open(vars_table, "r+") as infile:
        for line in infile:
            if header:
                print "Dif\t"+line
                header = False
            else:
                try:
                    linelist = line.split()
                    depth_threhold = 30

                    if ((float(linelist[11].split(",")[1]) + float(linelist[11].split(",")[0])) > depth_threhold and
                        (float(linelist[5].split(",")[1]) + float(linelist[5].split(",")[0])) > depth_threhold and
                        (float(linelist[8].split(",")[1]) + float(linelist[8].split(",")[0])) > depth_threhold and
                        (float(linelist[14].split(",")[1]) + float(linelist[14].split(",")[0])) > depth_threhold):


                        ratio_3D7 = float(linelist[5].split(",")[0]) / (float(linelist[5].split(",")[1]) + float(linelist[5].split(",")[0]))
                        ratio_B11 = float(linelist[8].split(",")[0]) / (float(linelist[8].split(",")[1]) + float(linelist[8].split(",")[0]))
                        ratio_E5HA = float(linelist[11].split(",")[0]) / (float(linelist[11].split(",")[1]) + float(linelist[11].split(",")[0]))
                        ratio_NF54 = float(linelist[14].split(",")[0]) / (float(linelist[14].split(",")[1]) + float(linelist[14].split(",")[0]))

                        #print ratio_3D7, ratio_B11, ratio_E5HA, ratio_NF54
                        contrasts = [abs(ratio_E5HA-ratio_NF54), abs(ratio_3D7-ratio_B11), abs(ratio_3D7-ratio_E5HA), abs(ratio_3D7-ratio_NF54), abs(ratio_3D7-ratio_B11),
                                    abs(ratio_B11-ratio_E5HA), abs(ratio_B11-ratio_NF54)]

                        difs = ["E5HA_NF54", "3D7_B11", "3D7_E5HA", "3D7_NF54", "3D7_B11", "B11_E5HA", "B11_NF54"]

                        # print contrasts
                        # print "\n"

                        if any(i > 0.8 for i in contrasts):
                            mask = [x > 0.8 for x in contrasts]
                            print "Differential variant:    ", list(compress(difs, mask)), "\t", line,
                            #print contrasts


                except:
                    pass
                    #print "Bad line: ", line


if __name__ == "__main__":
	filenames = sys.argv[1:]
	for element in filenames:
		parse_variants(element)
