#!/usr/bin/env python

'''
Reformat facets cn files to a bedtools-compatible format (chr start stop)
Correct chromosome names if they have been altered by facets.

Franziska Singer, January 2019
'''

import sys
import numpy as np
import os
import re
import argparse
import math

'''
function definitions
'''

def checkChromNames():
	inFasta = open(args.refFile,'r')
	refLine = inFasta.readline().strip()
	inFasta.close()

	if not refLine.startswith(">"):
		print("Error. Reference fasta file is not properly formatted!")
		return ""

	if refLine[1:].startswith("chr"):  # Note: by default facets only check for chr in the sample and removes this string.
		return "chr"
	else:
		return ""

def getColumnIndex(firstInputLine):
	firstLineSplit = firstInputLine.strip().split('\t')
	index_chrom = -1
	index_start = -1
	index_stop = -1
	
	for pos in range(0,len(firstLineSplit)):
		if args.colName_chrom == firstLineSplit[pos]:
			index_chrom = pos
			continue
		if args.colName_start == firstLineSplit[pos]:
			index_start = pos
			continue
		if args.colName_stop == firstLineSplit[pos]:
			index_stop = pos
			continue
	
	if (index_chrom == -1) or (index_start == -1) or (index_stop == -1):
		print("Error! Did not find all columns matching the input parameteres in input header: \n%s." %(firstInputLine))
		sys.exit(1)
	
	return (index_chrom, index_start, index_stop)


'''
main body
'''

parser = argparse.ArgumentParser(description='Make facets .cn format bedtools compatible to allow CNV annotation. Further, check and correct if chromosome names have been altered by facets.')
parser.add_argument('--inFile', dest='inFile', required=True, help='Input table with copy number calls, tab separated.')
parser.add_argument('--refFile', dest='refFile', required=True, help='Fasta file with the reference sequence used for mapping.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the reformatted output file.')
parser.add_argument('--colName_chrom', dest='colName_chrom', required=True, help='Column name of column containing the chromosome names.')
parser.add_argument('--colName_start', dest='colName_start', required=True, help='Column name of column containing the CNV start position.')
parser.add_argument('--colName_stop', dest='colName_stop', required=True, help='Column name of column containing the CNV stop position.')
#parser.add_argument('--colName_minorCopy', dest='colName_minorCopy', required=True, help='Column name of column containing the minor alelle copy state.')

args = parser.parse_args()

# first check whether the reference fasta file contains chromosome names with "chr"

chrom_chr = checkChromNames() # empty string if nothing needs to be changed

infile = open(args.inFile,'r')
firstInputLine = infile.readline().strip()

# get comlumn indices

(index_chrom, index_start, index_stop) = getColumnIndex(firstInputLine)

outfile = open(args.outFile,'w')

# Necessary for output: 
# 1) header begins with '#'
# 2) columns 1,2,3: chrom, start, stop -> reorder features

#example of original input:
#chrom   seg     num.mark        nhet    cnlr.median     mafR    segclust        cnlr.median.clust       mafR.clust      start   end     cf.em   tcn.em  lcn.em  clonal.cluster

firstLineSplit = firstInputLine.strip().split('\t')
outfileHeader = "#%s\t%s\t%s" %(firstLineSplit[index_chrom],firstLineSplit[index_start],firstLineSplit[index_stop])

for pos in range(0,len(firstLineSplit)):
	if (pos == index_chrom) or (pos == index_start) or (pos == index_stop):
		continue
	outfileHeader = "%s\t%s" %(outfileHeader,firstLineSplit[pos])

outfile.write(outfileHeader + "\n")

for line in infile:
	lineSplit = line.strip().split('\t')
	chromName = lineSplit[index_chrom]
	if lineSplit[index_chrom] == "23":
		chromName = "X"

	lineTemp = "%s%s\t%s\t%s" %(chrom_chr,chromName,lineSplit[index_start],lineSplit[index_stop])

	for posi in range(0,len(lineSplit)):
		if (posi == index_chrom) or (posi == index_start) or (posi == index_stop):
			continue

		lineTemp = "%s\t%s" %(lineTemp,lineSplit[posi])
	
	outfile.write(lineTemp + "\n")

infile.close()
outfile.close()

