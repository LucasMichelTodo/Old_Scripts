#!/usr/bin/env python

import sys
import subprocess

with open("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Tcd4_predictions/tcd4_doubleHLA_predictions.txt", "r+") as infile:
  for line in infile:
    cmd = line.strip()
    subprocess.call(cmd, shell=True)
