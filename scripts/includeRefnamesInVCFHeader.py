#!/usr/bin/env python

'''
Integrate the reference contigs from a separate file to given vcf file (e.g. required for varscan or freebayes output).
Franziska Singer, August 2015
'''

import sys

if len(sys.argv) <= 1:
	print "Integrate the reference contigs from a separate file to given vcf file (e.g. required for varscan or freebayes output)."
	print "Usage: python includeRefnamesInVCFHeader.py [vcfFile] [referenceNamesList] [outFile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Integrate the reference contigs from a separate file to given vcf file (e.g. required for varscan or freebayes output)."
	print "Usage: python includeRefnamesInVCFHeader.py [vcfFile] [referenceNamesList] [outFile]"
	sys.exit(1)

vcfFile = sys.argv[1]
refNamesFile = sys.argv[2]
outFile = sys.argv[3]

refNameLines = []    # contains the reference contig names

infileRef = open(refNamesFile,'r')

for lineref in infileRef:
	refNameLines.append(lineref.strip())

infileRef.close()

outfile = open(outFile,'w')

infile = open(vcfFile,'r')

for line in infile:
	if line.startswith("##"):
		outfile.write(line)
		continue
	elif line.startswith("#"):
		for refNameTemp in refNameLines:
			outfile.write(refNameTemp + "\n")
		outfile.write(line)
		continue
	outfile.write(line)
	
	
infile.close()
outfile.close()
