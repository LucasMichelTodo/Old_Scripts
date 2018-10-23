#!/usr/bin/env python


### First import all probe annotations and convert into a dictionary:

probe_dict = {}

with open("/home/lucas/ISGlobal/Arrays/good_probes_final.csv", "r+") as good_probes:
    header = True
    for line in good_probes:
        if header == True:
            cols = line.strip().split(",")
            header = False
        else:
            probe = line.strip().split(",")
            for i in range(len(cols)):
                if i == 0:
                    probe_dict[probe[0]] = {}
                else:
                    probe_dict[probe[0]][cols[i]] = probe[i]

            probe_dict[probe[0]]["status"] = "keep"


with open("/home/lucas/ISGlobal/Arrays/excepcions_probes_final.csv", "r+") as exception_probes:
    header = True
    for line in exception_probes:
        if header == True:
            cols = line.strip().split(",")
            header = False
        else:
            probe = line.strip().split(",")
            for i in range(len(cols)):
                if i == 0:
                    probe_dict[probe[0]] = {}
                else:
                    probe_dict[probe[0]][cols[i]] = probe[i]

            probe_dict[probe[0]]["status"] = "keep"


with open("/home/lucas/ISGlobal/Arrays/bad_probes_final.csv", "r+") as bad_probes:
    header = True
    for line in bad_probes:
        if header == True:
            cols = line.strip().split(",")
            header = False
        else:
            probe = line.strip().split(",")
            for i in range(len(cols)):
                if i == 0:
                    probe_dict[probe[0]] = {}
                else:
                    probe_dict[probe[0]][cols[i]] = probe[i]

            probe_dict[probe[0]]["status"] = "drop"
            probe_dict[probe[0]]["New_Target"] = "NA"
            probe_dict[probe[0]]["New_Annot"] = "NA"
            probe_dict[probe[0]]["Seq"] = "NA"



## Import array description and complete with info stored in the probe dictionary:

with open("/home/lucas/ISGlobal/Arrays/Array_Annotation/array_description_nou.csv", "r+") as array:
    header = True
    for line in array:
        if header == True:
            cols = line.strip().split(",")
            header = False
            print "\t".join(cols).replace("\"", "")+"\t"+\
            'Score_1'+"\t"+'Hit_1'+"\t"+'Score_2'+"\t"+'Hit_2'+"\t"+"New_Target"+"\t"+"New_Annot"+"\t"+"Status"+"\t"+'Seq'+"\t"+'Kasfak'+"\t"+'Target_Kasfak_New_ID'+"\t"+'Transcripts'
        else:
            probe = line.strip().split(",")
            probe_name = probe[2].replace("\"", "")
            if probe_name in probe_dict.keys():
                print "\t".join(probe)+"\t"+probe_dict[probe_name]["Score_1"]+"\t"+probe_dict[probe_name]["Hit_1"]+"\t"+probe_dict[probe_name]["Score_2"]+"\t"+probe_dict[probe_name]["Hit_2"]+\
                "\t"+probe_dict[probe_name]["New_Target"]+"\t"+probe_dict[probe_name]["New_Annot"]+"\t"+probe_dict[probe_name]["status"]+"\t"+probe_dict[probe_name]["Seq"]+\
                "\t"+probe_dict[probe_name]["Kasfak"]+"\t"+probe_dict[probe_name]["Target_Kasfak_New_ID"]+"\t"+probe_dict[probe_name]["Transcripts"]
            else:
                print "\t".join(probe)+"\t"+"\t".join(["NA"]*11)
