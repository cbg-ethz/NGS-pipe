#!/usr/bin/env python

'''
Filter BICseq2 results (after genotyping.pl)
Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re


if "-h" in sys.argv[1]:
	print("Filter BICse2 results and reformat to .bed file for bedtools annotation.")
	print("Usage: python filterBICseq2.py [infileBICseq2] [outfile] [pvalueThreshold]")
	sys.exit(1)

infileBICseq2 = sys.argv[1]
outfileName = sys.argv[2]
pvalueThreshold = float(sys.argv[3])

infile = open(infileBICseq2,'r')
outfile = open(outfileName,'w')

outfile.write("#Chromosome\tstart\tend\tlog2.copyRatio\tpvalue\testimated copy number\n")
allCNVs = 0
filteredCNVs = 0
infile.readline() # skip first line

for line in infile:   # format: chrom	start	end	binNum	tumor	tumor_expect	normal	normal_expect	log2.copyRatio	log2.TumorExpectRatio	pvalue	pvalue.TumorVsExpected
	lineArr = line.strip().split("\t")
				
	pvalue = float(lineArr[10])
		
	allCNVs += 1
	if pvalue > pvalueThreshold:
		filteredCNVs += 1
		continue
	
	#copyNum = (float(lineArr[4])/float(lineArr[6])) * 2.0  # does not include normalization
	copyNum = 2.0 * (2**float(lineArr[8]))
	outfile.write("%s\t%s\t%s\t%s\t%s\t%.1f\n"%(lineArr[0],lineArr[1],lineArr[2],lineArr[8],lineArr[10],copyNum))

infile.close()
outfile.close()

print("Filtered %s of %s CNVs with threshold %s." %(filteredCNVs,allCNVs,pvalueThreshold))
