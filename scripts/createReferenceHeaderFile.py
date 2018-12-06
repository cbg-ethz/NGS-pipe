#!/usr/bin/env python

'''
Given a sam file header, create a txt file containing contig names and lengths necessary to update the vcf header (required e.g. for VarScan2 to enable GATK variantCombine)
Franziska Singer, February 2018
'''

import sys
import numpy as np
import os
import re


if (len(sys.argv) <= 1) or "-h" in sys.argv[1]:
    print("Create reference names vcf header file.")
    print("Usage: python createReferenceHeaderFile.py [samHeader] [outfile]")
    sys.exit(1)

samHeaderFile = sys.argv[1]
outfileName = sys.argv[2]
infile = open(samHeaderFile,'r')
outfile = open(outfileName,'w')

for line in infile:
    if not line.startswith("@SQ"):
        continue
    lineSplit = line.strip().split("\t")
    chrName = lineSplit[1].split(":")[1]
    length = lineSplit[2].split(":")[1]

    outfile.write("##contig=<ID=" + chrName + ",length=" + length + ">\n")

infile.close()
outfile.close()
