#!/usr/bin/env python

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/B_predictions/New_Methods/cbtope_predictions.txt") as filein:
  for line in filein:
      try:
          if (line.split()[3] != "Non-Epitope") and (line.split()[3] != "Score"):
              print line.split()[3]
      except:
          pass
