#!/usr/bin/env python

# Import packages 
import pandas

gff = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", sep="\t", header=None, index_col=None, skiprows=18)

bed = pandas.read_csv(filepath_or_buffer= "/home/lucas/ISGlobal/Chip_Seq/DATA/Substractions/Macs_diff/diff_A7_E5_narrowpeaks_c3.0_cond1.bed", sep="\t", header=None, index_col=None, skiprows=1)

print bed[0:5]

#for line in bed 

