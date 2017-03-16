import glob,sys

rootdir = sys.argv[1]
mode = sys.argv[2]#'rna_xeno' #rna|wes|wgs|rna_xeno
if mode == 'rna_xeno':
    fastqdir = rootdir + 'fastq/'
    cliptrimdir = rootdir + 'cliptrim1/'
    alndir = rootdir + 'alignHumanOnlyOut3/'
    fcntdir = rootdir + 'featurecounts4/'
    print 'SAMPLE\tFASTQ\tCLIPTRIM\tCLIPTRIM_PERCENTAGE\tALIGNED\tALIGNED_PERCENTAGE\tREADS_ON_GENES\tREADS_ON_GENES_PERCENTAGE\tGENES'

if mode == 'rna':
    fastqdir = rootdir + 'fastq/'
    cliptrimdir = rootdir + 'cliptrim1/'
    alndir = rootdir + 'align2/'
    fcntdir = rootdir + 'featurecounts3/'
    print 'SAMPLE\tFASTQ\tCLIPTRIM\tCLIPTRIM_PERCENTAGE\tALIGNED\tALIGNED_PERCENTAGE\tREADS_ON_GENES\tREADS_ON_GENES_PERCENTAGE\tGENES'
if mode =='wgs':
    fastqdir = rootdir + 'fastq/'
    cliptrimdir = rootdir + 'cliptrim/'
    alndir = rootdir + 'merged/'
    nodupdir = rootdir + 'removedPcrDuplicates/'
    print 'SAMPLE\tFASTQ\tCLIPTRIM\tCLIPTRIM_PERCENTAGE\tALIGNED\tALIGNED_PERCENTAGE\tNO_DUPLICATES\tNO_DUPLICATES_PERCENTAGE'

if mode == 'wes':
    fastqdir = rootdir + 'fastq/'
    cliptrimdir = rootdir + 'cliptrim/'
    alndir = rootdir + 'merged/'
    nodupdir = rootdir + 'removedPcrDuplicates/'
    recaldir = rootdir + 'recalibratedBases/'
    print 'SAMPLE\tFASTQ\tCLIPTRIM\tCLIPTRIM_PERCENTAGE\tALIGNED\tALIGNED_PERCENTAGE\tNO_DUPLICATES\tNO_DUPLICATES_PERCENTAGE\tRECALIBRATED_BASES\tRECALIBRATED_BASES_PERCENTAGE'

def getReadCountFromCountFiles(counts):
    cnt = 0
    for count in counts:
        infile = open(count, 'rU')
        cnt = cnt + int(infile.readline().strip())
        infile.close()
    return cnt
def getReadCountFromFlagstatFile(flagstat):
    infile = open(flagstat, 'rU')
    mapped = 0
    secondary = 0
    linecount = 0
    for line in infile:
        if 'in total' in line:
            alignments = int(line.split('+')[0].strip()) + int(line.split('+')[1].strip().split('in total')[0].strip())
        if 'secondary' in line:
            secondaryalignments = int(line.split('+')[0].strip()) + int(line.split('+')[1].strip().split('secondary')[0].strip())
        if 'mapped (' in line:
            mappedReads = int(line.split('+')[0].strip()) + int(line.split('+')[1].strip().split('mapped (')[0].strip())
    infile.close()
    return (alignments,secondaryalignments, mappedReads)
def ac(instr):
    instr = str(instr)
    out = [instr[0]]
    for i in range(1, len(instr)):
        if (len(instr) - i) % 3 == 0:
            out.append(',')
        out.append(instr[i])
    return ''.join(out)
def getFeatureCountsCount(fcnt):
    featurecounts = 0
    genes = 0
    featurecountsfile = open(fcnt, 'rU')
    for line in featurecountsfile:
        featurecounts = featurecounts + int(line.split()[1].strip())
        if int(line.split()[1].strip()) > 0:
            genes = genes + 1
    featurecountsfile.close()
    return (featurecounts, genes)

SAMPLENAMES = [directory.replace(cliptrimdir,'') for directory in glob.glob(cliptrimdir + '*')]

for sample in SAMPLENAMES:
    countfiles = glob.glob(fastqdir + '/'  + sample + '/*/*.count')
    readCountFastq = getReadCountFromCountFiles(countfiles)
    countfiles = glob.glob(cliptrimdir + '/'  + sample + '/*/*.count')
    readCountcliptrim = getReadCountFromCountFiles(countfiles)
    readCountcliptrimPerc = str("{0:.2f}".format(float(readCountcliptrim*100.0)/readCountFastq))
    if 'rna' in mode:
        alignTuple = getReadCountFromFlagstatFile(glob.glob(alndir + '/' + sample + '.flagstat')[0])
        alignedReads = alignTuple[2] - alignTuple[1]
        alignedReadsPerc = str("{0:.2f}".format(float(alignedReads*100.0)/readCountFastq))
        (featurecounts, genes) = getFeatureCountsCount(fcntdir + sample + '.fcnt2htseq.txt')
        featurecountsPerc = str("{0:.2f}".format(float(featurecounts*100.0)/readCountFastq))
        print sample + '\t' + ac(readCountFastq) + '\t' + ac(readCountcliptrim) + '\t' + readCountcliptrimPerc + '\t' + ac(alignedReads) + '\t' + alignedReadsPerc + '\t' + ac(featurecounts) + '\t' + featurecountsPerc + '\t' + ac(genes)
    else:
        alignTuple = getReadCountFromFlagstatFile(alndir + '/' + sample + '.bam.flagstat')
        alignedReads = alignTuple[2] - alignTuple[1]
        alignedReadsPerc = str("{0:.2f}".format(float(alignedReads*100.0)/readCountFastq))
        alignTuple = getReadCountFromFlagstatFile(nodupdir + '/' + sample + '.bam.flagstat')
        alignedReadsNoDupl = alignTuple[2] - alignTuple[1]
        alignedReadsNoDuplPerc = str("{0:.2f}".format(float(alignedReadsNoDupl*100.0)/readCountFastq))
        if 'wgs' in mode:
            print sample + '\t' + ac(readCountFastq) + '\t' + ac(readCountcliptrim) + '\t' + readCountcliptrimPerc + '\t' + ac(alignedReads) + '\t' + alignedReadsPerc + '\t' + ac(alignedReadsNoDupl) + '\t' + alignedReadsNoDuplPerc
        else:
            alignTuple = getReadCountFromFlagstatFile(recaldir + '/' + sample + '.bam.flagstat')
            alignedReadsRecal = alignTuple[2] - alignTuple[1]
            alignedReadsRecalPerc = str("{0:.2f}".format(float(alignedReadsRecal*100.0)/readCountFastq))
            print sample + '\t' + ac(readCountFastq) + '\t' + ac(readCountcliptrim) + '\t' + readCountcliptrimPerc + '\t' + ac(alignedReads) + '\t' + alignedReadsPerc + '\t' + ac(alignedReadsNoDupl) + '\t' + alignedReadsNoDuplPerc + '\t' + ac(alignedReadsRecal) + '\t' + alignedReadsRecalPerc 
