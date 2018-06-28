#!/usr/bin/env python

import sys
import numpy as np

b_factors = []
norm_b_factors = []
rsas = []

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/PDB_results/kmp11/Results/kmp11.B99990001.pdb", "r+") as pdb_in:
    for line in pdb_in:
        if line.startswith("ATOM"):
            if line.strip().split()[2] == "CA":
                if  line.strip().split()[9] == "C":
                    bfact = float(line.strip().split()[8][4:])
                else:
                    bfact = float(line.strip().split()[9])
                b_factors.append(bfact)

mean = np.mean(b_factors)
sd = np.std(b_factors)

for i in b_factors:
    norm_b_factors.append((i-mean)/sd)


with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/PDB_results/kmp11/Results/kmp11.rsa", "r+") as nacces_in:
    for line in nacces_in:
        if line.startswith("RES"):
            rsas.append(line.strip().split()[3])

print norm_b_factors
print rsas
print "....."
print len(norm_b_factors)
print len(rsas)

epitopes = []
for i in range(len(rsas)):
    if (norm_b_factors[i] >= 0.5) & (rsas[i] > 50):
        epitopes.append("E")
    else:
        epitopes.append("-")

print "".join(epitopes)
