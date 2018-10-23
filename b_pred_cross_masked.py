#!/usr/bin/env python

import sys

def cross_predictions_and_masked(pred_file, cons_file):

    pred_dict = {}
    with open(pred_file, "r+") as file1:
        for line in file1:
            if line.startswith(">"):
                prot = line.replace(">", "").strip()

            else:
                pred_dict[prot] = line.strip()


    #print pred_dict
    masked_dic = {}
    with open(cons_file, "r+") as file2:
        for line in file2:
            if line.startswith(">"):
                prot = line.replace(">", "").split("/")[0].strip()

            else:
                masked_dic[prot] = line.strip()

    for ag, maskseq in masked_dic.iteritems():

        i = 0
        epitope_seq = ""

        if ag in pred_dict.keys():

            print ">"+ag
            for nn in maskseq:
                if i < len(pred_dict[ag]):

                    if pred_dict[ag][i] != "-":
                        epitope_seq += nn
                        i += 1
                    else:
                        epitope_seq += "-"
                        i += 1

            print epitope_seq



if __name__ == "__main__":
    filename= sys.argv[1:]
    cross_predictions_and_masked(filename[0], filename[1])
