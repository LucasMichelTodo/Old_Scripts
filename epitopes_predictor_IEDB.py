#!/usr/bin/env python

import sys
import subprocess

def predict_epitopes(fasta_file):

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

    methods = ["ann", "comblib_sidney2008", "consensus", "IEDB_recommended", "netmhcpan", "smm", "smmpmbec", "pickpocket", "netmhccons", "netmhcstabpan"]

    for i in methods:
        cmd = "mkdir {}" .format(i)
        subprocess.call(cmd, shell=True)

    for i in hlas:
        for x in methods:
            cmd = "/home/lucas/Programs/mhc_i/src/predict_binding.py {} {} 9 {} > ./{}/{}_{}.txt" .format(x, i, fasta_file, x, i, x)
            print cmd
            try:
                subprocess.call(cmd, shell=True)
            except:
                pass




if __name__ == "__main__":
    filenames = sys.argv[1:]
    for file in filenames:
        predict_epitopes(file)
