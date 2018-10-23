#!/usr/bin/env python

import sys
import re

seq = []
seqs = {}
aas = {}
aa = []

## *********************************************************** Method 1 ****************************************************************************
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/predictions_antigens.txt", "r+") as infile:
    for line in infile:
        if line.startswith("input:"):
            #print seq
            try:
                seqs[antigen] = seq
                aas[antigen] = aa
            except:
                pass
            seq = []
            aa = []
            antigen = line.strip().split()[1]
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            pass
        else:
            seq.append(float(line.split("\t")[2]))
            aa.append(line.split("\t")[1])

    #print seq
    seqs[antigen] = seq
    aas[antigen] = aa

# for ag, sq in seqs.iteritems():
#     print ag, len(sq)
#
#
# print "-----------------------------------------------------------------------------------"

## *********************************************************** Method 2 ****************************************************************************
seq = []
seqs2 = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/Chou-Fasman/antigens_Chou-Fasman_predicitons.txt", "r+") as infile:
    for line in infile:
        #print line
        if line.startswith("input:"):
            #print seq
            try:
                seqs2[antigen] = seq
            except:
                pass
            seq = []
            antigen = line.strip().split()[1]
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            pass
        elif line.startswith("*"):
            pass
        else:
            seq.append(float(line.split("\t")[5]))
            #seq.append(line.split("\t")[1])

    seqs2[antigen] = seq

for ag, sq in seqs2.iteritems():
    seqs2[ag] = [0.0,0.0,0.0]+sq+[0.0,0.0,0.0]

## *********************************************************** Method 3 ****************************************************************************
# for ag, sq in seqs2.iteritems():
#     print ag, len(sq)
#
# print "------------------------------------------------------------------------------------------"
#

seq = []
seqs3 = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/All_methods/antigens_Emini_predicitons.txt", "r+") as infile:
    for line in infile:
        #print line
        if line.startswith("input:"):
            #print seq
            try:
                seqs3[antigen] = seq
            except:
                pass

            seq = []
            antigen = line.strip().split()[1]

            avoid_block = True
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            avoid_block = False
        elif line.startswith("*"):
            pass
        elif not avoid_block:
            seq.append(float(line.split("\t")[5]))
            #seq.append(line.split("\t")[1])

    seqs3[antigen] = seq

for ag, sq in seqs3.iteritems():
    seqs3[ag] = [0.0,0.0]+sq+[0.0,0.0,0.0]

## *********************************************************** Method 4 ****************************************************************************
seq = []
seqs4 = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/All_methods/antigens_Karplus-Schulz_predicitons.txt", "r+") as infile:
    for line in infile:
        #print line
        if line.startswith("input:"):
            #print seq
            try:
                seqs4[antigen] = seq
            except:
                pass
            seq = []
            antigen = line.strip().split()[1]
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            pass
        elif line.startswith("*"):
            pass
        else:
            seq.append(float(line.split("\t")[5]))
            #seq.append(line.split("\t")[1])

    seqs4[antigen] = seq

for ag, sq in seqs4.iteritems():
    seqs4[ag] = [0.0,0.0,0.0]+sq+[0.0,0.0,0.0,0.0]

## *********************************************************** Method 5 ****************************************************************************

seq = []
seqs5 = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/All_methods/antigens_Kolaskar-Tongaonkar_predicitons.txt", "r+") as infile:
    for line in infile:
        #print line
        if line.startswith("input:"):
            #print seq
            try:
                seqs5[antigen] = seq
            except:
                pass

            seq = []
            antigen = line.strip().split()[1]

            avoid_block = True
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            avoid_block = False
        elif line.startswith("*"):
            pass
        elif not avoid_block:
            seq.append(float(line.split("\t")[5]))
            #seq.append(line.split("\t")[1])

    seqs5[antigen] = seq

for ag, sq in seqs5.iteritems():
    seqs5[ag] = [0.0,0.0,0.0]+sq+[0.0,0.0,0.0]

## *********************************************************** Method 6 **************************************************************************** ##
seq = []
seqs6 = {}
with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/All_methods/antigens_Parker_predicitons.txt", "r+") as infile:
    for line in infile:
        #print line
        if line.startswith("input:"):
            #print seq
            try:
                seqs6[antigen] = seq
            except:
                pass
            seq = []
            antigen = line.strip().split()[1]
            #print ">{}" .format(antigen)
        elif line.startswith("Position"):
            pass
        elif line.startswith("*"):
            pass
        else:
            seq.append(float(line.split("\t")[5]))
            #seq.append(line.split("\t")[1])

    seqs6[antigen] = seq

for ag, sq in seqs6.iteritems():
    seqs6[ag] = [0.0,0.0,0.0]+sq+[0.0,0.0,0.0]

## *********************************************************** Put it all toghether ****************************************************************************

scores = {}
masked_seqs = []

for ag in aas.keys():
    scores[ag] = zip(seqs[ag], seqs2[ag], seqs3[ag], seqs4[ag], seqs5[ag], seqs6[ag])

    seq = ""
    for i in range(1,len(aas[ag])):
        if scores[ag][i][0] >= 0.6 and all(x > 0.0 for x in scores[ag][i][1:5]) and scores[ag][i][5] >= -99999:
            seq += aas[ag][i]
        else:
            seq += "-"

    masked_seqs.append((ag,seq))

for i in masked_seqs:
    print i[0]
    print i[1]

print "********************************************************************************************"


with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/epitopes.fasta", "r+") as infile2:
    for line in infile2:
        if line.startswith(">"):
            pass
        else:
            regex = re.compile(re.escape(str(line.strip())))
            for i in masked_seqs:
                result = re.findall(regex, i[1])
                if result:
                    print result, i[0]
