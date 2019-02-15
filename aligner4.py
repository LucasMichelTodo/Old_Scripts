#!/usr/bin/env python

import sys
import subprocess


def align(filenames):

    subprocess.call("mkdir Results", shell=True)

    tupple_list = []
    i = 0
    while i < len(filenames):
        tupple_list.append((filenames[i], filenames[i+1]))
        i += 2

    for tupple in tupple_list:

        cmd = "~/Programs/bowtie2-2.3.0-legacy/bowtie2 -p 4 -t --very-sensitive --local -5 4 -3 4 -I 50 -X 200 -x ../Reference_Genomes/Ind_NI_Revised -1 {} -2 {} > ./Results_edited_genomes/{}"\
            .format(tupple[0],
                    tupple[1],
                    tupple[0].replace("_read1_clean.fastq.gz", "NI_Revised.sam"))

        subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    filenames = sys.argv[1:]
    align(filenames)
