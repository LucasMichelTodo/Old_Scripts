#!/usr/bin/env python

from collections import defaultdict

input_dir = "/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Run_27_02_18/fastas/iedb_methods/"

methods = ["smmpmbec", "pickpocket", "netmhccons", "netmhcstabpan", "ann", "smm", "comblib_sidney2008", "netmhcpan"]
predictions = defaultdict(list)

hlas = [
"HLA-A*02:01",
"HLA-A*02:02",
"HLA-A*02:05",
"HLA-A*68:02",
"HLA-A*03:01",
"HLA-A*11:01",
"HLA-A*31:01",
"HLA-A*33:01",
"HLA-A*68:01",
"HLA-A*66:01",
"HLA-A*24:02",
"HLA-B*38:01",
"HLA-B*07:02",
"HLA-B*35:01",
"HLA-B*51:01",
"HLA-B*51:02",
"HLA-B*53:01",
"HLA-A*01:01",
"HLA-B*15:01"
]

for i in methods:
    for x in hlas:
        header = True
        try:
            with open(input_dir+i+"/"+x+"_"+i+".txt", "r+") as filein:
                for line in filein:
                    if header:
                        keys = line.strip().split("\t")
                        with open("_".join(["./Selections/", "selection", i, x, ".txt"]), "w+") as outfile:
                            outfile.write("\t".join(keys)+"\n")

                        header = False
                    else:
                        values = line.strip().split("\t")
                        results = dict(zip(keys, values))
                        if float(results["rank"]) <= 0.1:
                            for a in values:
                                with open("_".join(["./Selections/", "selection", i, x, ".txt"]), "a+") as outfile:
                                    outfile.write(a+"\t")
                            with open("_".join(["./Selections/", "selection", i, x, ".txt"]), "a+") as outfile:
                                outfile.write("\n")
        except:
            print "Method {} has no HLA: {}!" .format(i,x)
# for x in hlas:
#     header = True
#     with open(input_dir+"IEDB_recommended"+"/"+x+"_"+"IEDB_recommended"+".txt", "r+") as filein:
#         for line in filein:
#             if header:
#                 keys = line.strip().split("\t")
#                 iedb_recommended = ["ann", "smm", "comblib_sidney2008", "netmhcpan"]
#                 for m in iedb_recommended:
#                     with open("_".join(["./Selections/", "selection", m, x, ".txt"]), "w+") as outfile:
#                         outfile.write("\t".join(keys[0:6])+"\t"+"ic50"+"\t"+"rank"+"\n")
#
#                 header = False
#             else:
#                 values = line.strip().split("\t")
#                 results = dict(zip(keys, values))
#                 for m in iedb_recommended:
#                     rank = m+"_rank"
#                     if results[rank] != "-":
#                         if float(results[rank]) <= 0.1:
#                             print line
                #
                # if float(results["ann_rank"]) <= 0.1:
                #     for a in values:
                #         with open("_".join(["./Selections/", "selection", "IEDB_recommended", x, ".txt"]), "a+") as outfile:
                #             outfile.write(a+"\t")
                #     with open("_".join(["./Selections/", "selection", "IEDB_recommended", x, ".txt"]), "a+") as outfile:
                #        outfile.write("\n")



                    #predictions[i].append(results)


    # for key, value in predictions.iteritems():
    #     print key
        # for x in value:
        #     if float(x["rank"]) < 2:
        #         for a, b in x.iteritems():
        #             print "{}: {}" .format(a, b)
