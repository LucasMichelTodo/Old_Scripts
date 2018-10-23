#!/usr/bin/env python

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Edit this line with output of "tcd4_predictions_parse.py" as parsed_sorted_cores and prediction files:

parsed_sorted_cores = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/IEDB_recommended/Tcd4_predictions/all_tcd4predictions_sorted.csv"
original_predictions = "/home/lucas/ISGlobal/Brucei/aat_vaccine/Run_060818/IEDB_recommended/Tcd4_predictions/{}_predictions.txt"


with open(parsed_sorted_cores, "r+") as infile:
    for line in infile:
        linelist = line.split("\t")

        if linelist[1] > 0:
            core = linelist[0]
            files = linelist[2].replace("[", "").replace("./", "").replace("'", "").replace("]", "").replace(",", "").split()

            seqs = []

            for x in files:
                with open(original_predictions .format(x), "r+") as readfile:
                    firstline = True
                    for line in readfile:
                        if firstline:
                            #print line.split()[5]
                            firstline = False
                        else:
                            linelist2 = line.split()
                            if linelist2[5] == "COMB.LIB.-SMM-NN":
                                if linelist2[7] == core:
                                    seqs.append(linelist2[4])
                            elif linelist2[5] == "SMM-NN-Sturniolo":
                                if linelist2[10] == core:
                                    seqs.append(linelist2[4])
                            elif linelist2[5] == "NetMHCIIpan":
                                if linelist2[16] == core:
                                    seqs.append(linelist2[4])
                            elif linelist2[5] == "SMM":
                                if linelist2[10] == core:
                                    seqs.append(linelist2[4])

            # From a list of strings ("ASTSV", "JAHGAS", ...) create a list of SeqRecord objects:
            fseqs = []
            y = 1
            for i in seqs:
                fseqs.append(SeqRecord(Seq(i), id = "Seq_{}" .format(y)))
                y += 1

            # If more than one sequence is present for one core, proceed to align them and obtain consensus:
            if len(fseqs) > 1:

                # Create a handle for a virtual fasta file
                handle = StringIO()
                # Write list of SeqObjects to the hadle in fasta format
                SeqIO.write(fseqs, handle, "fasta")
                # Read handle to data
                data = handle.getvalue()

                # Create a call to muscle with clustalw output format as only parameter
                muscle_cline = MuscleCommandline("muscle3.8.31_i86linux64", clwstrict=True)

                # Make the call with data as input and capture output and stderr
                stdout, stderr = muscle_cline(stdin=data)

                # Read stdout as string to an alignment object
                align = AlignIO.read(StringIO(stdout), "clustal")

                # Use summary and dumb consensus to retrieve whole consensus sequence
                summary_align = AlignInfo.SummaryInfo(align)
                final_seq = summary_align.dumb_consensus(0.01)

                print core, final_seq

            else:
                print fseqs[0].seq, fseqs[0].seq
