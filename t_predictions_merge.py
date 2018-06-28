#!/usr/bin/env python

import sys
from collections import Counter

def merge_predictions(predictions_list):

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

    for hla in hlas:
        dictionary = {}
        for prediction_file in predictions_list:
            if hla in prediction_file:
                with open(prediction_file, "r+") as infile:
                    header = True
                    for line in infile:
                        if header:
                            header = False
                        else:
                            method = str(prediction_file).split("_")[2]
                            epitope = line.split()[5]
                            dictionary.setdefault(method, []).append(epitope)

        #print "*****************************************{}**************************************************************" .format(hla)
        setlist = []
        for x in dictionary.values():
            setlist.append(set(x))
        u = set.intersection(*setlist)
        #print u
        flat_list = [item for sublist in setlist for item in sublist]
        counts = Counter(flat_list)
        for epi in sorted(counts, key=counts.get, reverse=True):
            if counts[epi] > 4:
                print epi#, counts[epi]





    # epitopes = {}
    # for prediction_file in predictions_list:
    #     with open(prediction_file, "r+") as infile:
    #         method = str(prediction_file).split("_")[2]
    #         for line in infile:
    #             epitopes.setdefault(line.split()[5],[]).append(method)
    #
    # for epitope, methods in epitopes.iteritems():
    #     if len(set(methods)) > 1:
    #         print "Epitope: {} >>> Methods {}" .format(epitope, set(methods))



if __name__ == "__main__":
    filenames = sys.argv[1:]
    merge_predictions(filenames)
