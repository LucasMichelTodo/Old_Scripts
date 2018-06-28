#!/usr/bin/env python
import sys
import subprocess

with open("/home/lucas/ISGlobal/Arrays/084561_D_SequenceList_20161220.txt", "r+") as infile:
  for line in infile:
      with open("/home/lucas/ISGlobal/Arrays/Probes_fastas/"+line.strip().split()[0].replace("/", ":")+".fasta", "w+") as outfile:
          outfile.write(">"+line.strip().split()[0]+"\n"+line.strip().split()[1])


blastn -query probes.fasta -subject ../Gen_Referencies/PlasmoDB-30_Pfalciparum3D7_AnnotatedTranscripts.fasta -outfmt 6 > balst1.txt
      # print ">"+line.strip().split()[0]
      # print line.strip().split()[1]
