#!/usr/bin/env python

import sys
import subprocess

#methods = ["Chou-Fasman", "Emini", "Karplus-Schulz", "Kolaskar-Tongaonkar", "Parker"]
methods = ["Bepipred2"]

def multi_method_b_predict(infile):
    for method in methods:
        cmd = "mkdir {}" .format(method)
        subprocess.call(cmd, shell=True)
        cmd = "mkdir {}" .format("./"+method+"/Plots")
        subprocess.call(cmd, shell=True)
        cmd = "python ~/Programs/bcell_standalone/predict_antibody_epitope.py -m {} -f {} > {} --plot {}"\
        .format(method, infile, "./"+method+"/"+infile.replace(".fasta", "_"+method+"_predicitons.txt"), "./"+method+"/Plots")
        subprocess.call(cmd, shell=True)


if __name__ == "__main__":
	filenames= sys.argv[1:]
	for file in filenames:
		multi_method_b_predict(file)
