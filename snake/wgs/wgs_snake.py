import os, glob, sys

# functionality like getSampleNames
include: "../common/misc/misc_snake.py"

# Check if the uses specified the proper input and output directories
if not 'FASTQDIR' in globals():
    print('You have to specify the root directory of the fastq files!')
    sys.exit(1)
if not 'OUTDIR' in globals():
    print('You have to specify the root directory where the results will be generated!')
    sys.exit(1)
if not 'TMPDIR' in globals():
    print('You have to specify the root directory where temporary files will be stored!')
    sys.exit(1)

# This is the default order in which the programs are executed
# If the user specified a different order the user specified version is chosen.
if not 'CLIPTRIMIN' in globals():
    CLIPTRIMIN = FASTQDIR
if not 'CLIPTRIMOUT' in globals():
    CLIPTRIMOUT = OUTDIR + 'cliptrim/'
if not 'BWAIN' in globals():
    BWAIN = CLIPTRIMOUT
if not 'BWAOUT' in globals():
    BWAOUT = OUTDIR + 'bwa/'
if not 'FIXMATEANDSORTIN' in globals():
    FIXMATEANDSORTIN = BWAOUT
if not 'FIXMATEANDSORTOUT' in globals():
    FIXMATEANDSORTOUT = OUTDIR + 'fix_sorted/'
if not 'MERGEBAMSIN' in globals():
    MERGEBAMSIN = FIXMATEANDSORTOUT
if not 'MERGEBAMSOUT' in globals():
    MERGEBAMSOUT = OUTDIR + 'merged/'
if not 'NOSECONDARYALNIN' in globals():
    NOSECONDARYALNIN = MERGEBAMSOUT
if not 'NOSECONDARYALNOUT' in globals():
    NOSECONDARYALNOUT = OUTDIR + 'noSecondaryAln/'
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = NOSECONDARYALNOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markedDuplicates/'
if not 'REMOVEPCRDUBLICATESIN' in globals():
    REMOVEPCRDUBLICATESIN = MARKPCRDUBLICATESOUT
if not 'REMOVEPCRDUBLICATESOUT' in globals():
    REMOVEPCRDUBLICATESOUT = OUTDIR + 'removedPcrDuplicates/'
if not 'REALIGNINDELSIN' in globals():
    REALIGNINDELSIN = REMOVEPCRDUBLICATESOUT
if not 'REALIGNINDELSOUT' in globals():
    REALIGNINDELSOUT = OUTDIR + 'realignedIndels/'
if not 'BASERECALIBRATIONIN' in globals():
    BASERECALIBRATIONIN = REALIGNINDELSOUT
if not 'BASERECALIBRATIONOUT' in globals():
    BASERECALIBRATIONOUT = OUTDIR + 'recalibratedBases/'
if not 'BICSEQ2IN' in globals():
    BICSEQ2IN = REMOVEPCRDUBLICATESOUT
if not 'BICSEQ2OUT' in globals():
    BICSEQ2OUT = OUTDIR + 'bicseq2/'
if not 'MPILEUPIN' in globals():
    MPILEUPIN = REMOVEPCRDUBLICATESOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'
if not 'VARSCANCNVIN' in globals():
    VARSCANCNVIN = MPILEUPOUT
if not 'VARSCANCNVOUT' in globals():
    VARSCANCNVOUT = OUTDIR + 'varscan_cnv/'

# Definition of some constantly used lists
SINGLEFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
PAIREDFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
PAIREDFASTQFILESORPAHNS = [file.replace(FASTQDIR, '').replace('PAIREDEND','PAIREDEND/ORPHAN').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
PAIREDFASTQFILESWITHOUTR = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
PAIREDUNPAIRED = ['_PAIRED','_UNPAIRED']
MATEPAIR = ['_R1', '_R2']
SAMPLENAMES=getSampleNames()

# Include the rules
include: "../common/clip_trim/clip_trim_snake.py"
include: "../common/align/align_snake.py"
include: "../common/bam_mod/bam_mod_snake.py"
include: "../common/stats/stats_snake.py"
include: "../common/copy_number/copy_number_snake.py"

# This rule defines which files should be created
rule wgs:
    input:    
        expand(MERGEBAMSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(BICSEQ2OUT + '{tumorNormalMatching}.filtered.annotated.hg19_multianno.txt', tumorNormalMatching = getNormalTumorFiles())
    output:
        OUTDIR + 'complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log',
        scratch = 1000,
        mem = 1000,
        time = 10
    benchmark:
        OUTDIR + 'complete.txt.benchmark'
    shell:
        'date > {output}'

