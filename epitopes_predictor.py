#!/usr/bin/env python

import sys
import subprocess

def predict_epitopes(fasta_file):

    hlas = [
    "HLA-A0201",
    "HLA-A0202",
    "HLA-A0203",
    "HLA-A0205",
    "HLA-A0206",
    "HLA-A0207",
    "HLA-A6802",
    "HLA-A0301",
    "HLA-A1101",
    "HLA-A3101",
    "HLA-A3301",
    "HLA-A6801",
    "HLA-A6601",
    "HLA-A2402",
    "HLA-B3801",
    "HLA-B0702",
    "HLA-B3501",
    "HLA-B5101",
    "HLA-B5102",
    "HLA-B5301",
    "HLA-B5401",
    "HLA-A0101",
    "HLA-B1501",
    "HLA-B1502",
    ]

    for i in hlas:
        prefix = i.replace("-", "_")
        cmd = "netMHC -a {} -f {} > {}_binders.txt" .format(i, fasta_file, prefix)
        print cmd
        try:
            subprocess.call(cmd, shell=True)
        except:
            pass
        cmd = "awk \'$16 == \"SB\" {{ print $0 }}\' {}_binders.txt > {}_binders_SB.txt" .format(prefix, prefix)
        print cmd
        try:
            subprocess.call(cmd, shell=True)
        except:
            pass
        cmd = "awk '!/X/' {}_binders_SB.txt > {}_binders_SB_filetered.txt" .format(prefix, prefix)
        print cmd
        try:
            subprocess.call(cmd, shell=True)
        except:
            pass

    #awk '$16 == "SB" { print $0 }' netMHC_binders.txt > netMHC_binders_SB.txt
    #awk '!/X/' netMHC_binders_SB.txt > netMHC_binders_SB_filtered.txt

if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        predict_epitopes(file)
