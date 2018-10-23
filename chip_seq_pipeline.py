#!/usr/bin/env python

import sys
import os
import subprocess

def pipeline(filelist):

    str_filelist = " ".join(filelist)
    cmd = "aligner.py {}" .format(str_filelist)
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    filenames = sys.argv[1:]
    pipeline(filenames)
