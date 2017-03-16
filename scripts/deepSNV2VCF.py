
'''
Convert the output from deepSNV into vcf file format.
'''

import sys
import numpy as np
import os
import math
from Bio import SeqIO

if len(sys.argv) <= 1:
    print("\nConvert the output from deepSNV into vcf file format.\n")
    print("Usage: python deepSNV2VCF.py [deepSNV_variantList.txt] [reference_genome.fasta] [out.vcf]\n")
    print("Optional arguments, which are off by default:")
    print("\t--max-adj-pvalue <float>")
    print("\t--no-indels")
    print("\t--no-LOH")
    sys.exit(1)
if "-h" in sys.argv[1]:
    print("\nConvert the output from deepSNV into vcf file format.\n")
    print("Usage: python deepSNV2VCF.py [deepSNV_variantList.txt] [reference_genome.fasta] [out.vcf]\n")
    print("Optional arguments, which are off by default:")
    print("\t--max-adj-pvalue <float>")
    print("\t--no-indels")
    print("\t--no-LOH")
    sys.exit(1)


########
#functions
########
# this function gets the sequence neighborhood of a genomic position
def getNeighbors(chromosome, pos, records, num_neighbors):
        startPosBed=pos-(num_neighbors+1)   # 0-based
        endPosBed=pos+num_neighbors         # 1-based

        for record in records:
                if record.id == chromosome:
                        seqArr=list(record.seq[startPosBed:endPosBed])

        #print(seqArr)
        return seqArr



deepSNVList = sys.argv[1]
genomeFasta = sys.argv[2]
outName = sys.argv[3]
maxPval=1.0
indelsIn=True
lohIn=True
if len(sys.argv)>3:
    for i in range(3,len(sys.argv)):
        if sys.argv[i]=="--max-adj-pvalue":
            maxPval=float(sys.argv[i+1])
        if sys.argv[i]=="--no-indels":
            indelsIn=False
        if sys.argv[i]=="--no-LOH":
            lohIn=False



handle = open(genomeFasta, "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

infile = open(deepSNVList,'r')
if os.path.exists(outName):
    print("Outfile already exists: %s")%(outName)
    sys.exit(1)
outfile = open(outName,'w')

allVariants = 0
filteredVariants = 0

# print header lines already
outfile.write("##fileformat=VCFv4.0\n")
outfile.write("##source=deepSNV\n")
outfile.write("##reference=%s\n" % genomeFasta)
outfile.write("##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in normal\">\n")
outfile.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Coverage in normal\">\n")
outfile.write("##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in tumor\">\n")
outfile.write("##INFO=<ID=CT,Number=1,Type=Integer,Description=\"Coverage in tumor\">\n")
outfile.write("##INFO=<ID=VF,Number=.,Type=Float,Description=\"Frequency of the SNV in tumor\">\n")
outfile.write("##INFO=<ID=DVF,Number=.,Type=Float,Description=\"Relative frequency of the SNV from deepSNV\">\n")
outfile.write("##INFO=<ID=DSF,Number=.,Type=Float,Description=\"Estimated variance of the frequency from deepSNV\">\n")
outfile.write("##INFO=<ID=RP,Number=.,Type=Float,Description=\"Raw p-value\">\n")
outfile.write("##INFO=<ID=AP,Number=.,Type=Float,Description=\"Adjusted p-value\">\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")



# loop over list with variants
for line in infile:
    lineSplit = line.strip().split("\t")
    if lineSplit[1]=="pos":
        continue # header line
    allVariants += 1

    chr=lineSplit[0]
    posRaw=float(lineSplit[1]) # should be int, but there was a case in deepSNV like: chrX    4.1e+07 T       G
    ref=lineSplit[2]
    var=lineSplit[3]
    raw_pVal=float(lineSplit[4])
    freq_var=float(lineSplit[5])
    sigma2_freq_var=float(lineSplit[6])
    n_tst_fw=int(lineSplit[7])
    cov_tst_fw=int(lineSplit[8])
    n_tst_bw=int(lineSplit[9])
    cov_tst_bw=int(lineSplit[10])
    n_ctrl_fw=int(lineSplit[11])
    cov_ctrl_fw=int(lineSplit[12])
    n_ctrl_bw=int(lineSplit[13])
    cov_ctrl_bw=int(lineSplit[14])
    adjusted_pVal=float(lineSplit[15])

    pos=int(posRaw)

    n_tst=n_tst_fw + n_tst_bw
    cov_tst=cov_tst_fw + cov_tst_bw
    n_ctrl=n_ctrl_fw + n_ctrl_bw
    cov_ctrl=cov_ctrl_fw + cov_ctrl_bw

    if raw_pVal > 0:
        qual = -10 * math.log(raw_pVal,10)
    else:
        qual = 255

    VAF_normal = (float(n_ctrl)/cov_ctrl)
    VAF_tumor = (float(n_tst)/cov_tst)
    VAF_diff = VAF_tumor - VAF_normal
    isLOH=False
    if VAF_diff < 0.0:
        isLOH=True


    if adjusted_pVal <= maxPval:
        if lohIn or (not lohIn and not isLOH): # if we allow LOH, or we don't want LOH, but the current variant is not an LOH event
            # now check whether it is an indel or not
            if ref=="-" or var=="-": # indel
                if indelsIn:
                    numNeighbors=1
                    # get reference base to the left of it
                    baseNeighborhood=getNeighbors(chr, pos, records,numNeighbors) # the returned array will be [X,ref,X] the one-base-neighborhood of the ref allele
                    pos -= 1
                    if ref=="-": # insertion
                        ref=baseNeighborhood[0]
                        var=baseNeighborhood[0]+var
                    elif var=="-": # deletion
                        ref=baseNeighborhood[0]+baseNeighborhood[1]
                        var=baseNeighborhood[0]
                    # now write it out
                    outfile.write(
                        '%s\t%d\t.\t%s\t%s\t%.2f\t.\tNV=%d;CN=%d;TV=%d;CT=%d;VF=%.5f;DVF=%.5f;DSF=%.5f;RP=%f;AP=%f\n' % (chr, pos, ref, var, qual, n_ctrl, cov_ctrl, n_tst, cov_tst, VAF_tumor, freq_var, sigma2_freq_var, raw_pVal, adjusted_pVal))
                else:
                    filteredVariants+=1
                    #print("Filter because is indel, and indelsIn=%s" % (indelsIn))
                    continue
            else: # no indel, just write it out
                outfile.write(
                    '%s\t%d\t.\t%s\t%s\t%.2f\t.\tNV=%d;CN=%d;TV=%d;CT=%d;VF=%.5f;DVF=%.5f;DSF=%.5f;RP=%f;AP=%f\n' % (chr, pos, ref, var, qual, n_ctrl, cov_ctrl, n_tst, cov_tst, VAF_tumor, freq_var, sigma2_freq_var, raw_pVal, adjusted_pVal))
        else:
            filteredVariants+=1
            #print("Filter because of LOH: lohIN=%s, but isLOH=%s" % (lohIn, isLOH))
            continue
    else:
        filteredVariants+=1
        #print("Filter because of adjusted pvalue: %f" % (adjusted_pVal))
        continue


infile.close()
outfile.close()
percentageFiltered=float((float(filteredVariants)/float(max(1,allVariants)))*100.0)
print("Converted variants from the %s file into vcf file format. Filtered %s of %s variants (%.2f%%).") %(os.path.basename(deepSNVList),filteredVariants,allVariants,percentageFiltered)
