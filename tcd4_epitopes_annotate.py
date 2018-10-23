#!/usr/bin/env python

import re

fasta_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/all_proteomes.fasta"
epitopes_file = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/IEDB_recommended/Tcd4_predictions/tcd4_cores_and_consensus.csv"

for line in epitopes_file:
    regex = line.split()[1].strip
