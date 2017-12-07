import os, glob, sys
from itertools import chain

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
if not 'TRIMMOMATICIN' in globals():
    TRIMMOMATICIN = FASTQDIR
if not 'TRIMMOMATICOUT' in globals():
    TRIMMOMATICOUT = OUTDIR + 'cliptrim/'
if not 'BWAIN' in globals():
    BWAIN = TRIMMOMATICOUT
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
if not 'FREEBAYESIN' in globals():
    FREEBAYESIN = REALIGNINDELSOUT
if not 'FREEBAYESOUT' in globals():
    FREEBAYESOUT = OUTDIR + 'freebayes/raw/'
if not 'MPILEUPIN' in globals():
    MPILEUPIN = BASERECALIBRATIONOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'
if not 'VARSCANSNPIN' in globals():
    VARSCANSNPIN = MPILEUPOUT
if not 'VARSCANSNPOUT' in globals():
    VARSCANSNPOUT = OUTDIR + 'variants/varscansnp/raw/'
if not 'VARSCANSOMATICIN' in globals():
    VARSCANSOMATICIN = MPILEUPOUT
if not 'VARSCANSOMATICOUT' in globals():
    VARSCANSOMATICOUT = OUTDIR + 'variants/varscan_somatic/raw/'
if not 'VARSCANUPDATEHEADERIN' in globals():
    VARSCANUPDATEHEADERIN = VARSCANSOMATICOUT
if not 'VARSCANUPDATEHEADEROUT' in globals():
    VARSCANUPDATEHEADEROUT = OUTDIR + 'variants/varscan_somatic/complete_raw/'
if not 'VARSCANCOMPLETEIN' in globals():
    VARSCANCOMPLETEIN = VARSCANUPDATEHEADEROUT
if not 'VARSCANCOMPLETEOUT' in globals():
    VARSCANCOMPLETEOUT = OUTDIR + 'variants/varscan_somatic/combined_raw/'
if not 'VARSCANSOMATICFILTERIN' in globals():
    VARSCANSOMATICFILTERIN = VARSCANCOMPLETEOUT
if not 'VARSCANSOMATICFILTEROUT' in globals():
    VARSCANSOMATICFILTEROUT = OUTDIR + 'variants/varscan_somatic/filtered/'
if not 'BCFTOOLSIN' in globals():
    BCFTOOLSIN = MPILEUPOUT
if not 'BCFTOOLSOUT' in globals():
    BCFTOOLSOUT = OUTDIR + 'bcftools/raw/'
if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = BASERECALIBRATIONOUT
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = OUTDIR + 'variants/GATK/raw/'
if not 'MUTECT2IN' in globals():
    MUTECT2IN = BASERECALIBRATIONOUT
if not 'MUTECT2OUT' in globals():
    MUTECT2OUT = OUTDIR + 'variants/mutect2/'
if not 'MUTECT1IN' in globals():
    MUTECT1IN = BASERECALIBRATIONOUT
if not 'MUTECT1OUT' in globals():
    MUTECT1OUT = OUTDIR + 'variants/mutect1/raw/'
if not 'MUTECT1FILTERIN' in globals():
    MUTECT1FILTERIN = MUTECT1OUT
if not 'MUTECT1FILTEROUT' in globals():
    MUTECT1FILTEROUT = OUTDIR + 'variants/mutect1/filtered/'
if not 'JOINTSNVMIX2_075_IN' in globals():
    JOINTSNVMIX2_075_IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2_075_OUT' in globals():
    JOINTSNVMIX2_075_OUT = OUTDIR + 'variants/jointSNVMix2_075/raw/'
if not 'JOINTSNVMIX2IN' in globals():
    JOINTSNVMIX2IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2OUT' in globals():
    JOINTSNVMIX2OUT = OUTDIR + 'variants/jointSNVMix2/raw/'
if not 'SOMATICSNIPERIN' in globals():
    SOMATICSNIPERIN = BASERECALIBRATIONOUT
if not 'SOMATICSNIPEROUT' in globals():
    SOMATICSNIPEROUT = OUTDIR + 'variants/somaticSniper/raw/'
if not 'VARDICTIN' in globals():
    VARDICTIN = BASERECALIBRATIONOUT
if not 'VARDICTOUT' in globals():
    VARDICTOUT = OUTDIR + 'variants/varDict/raw/'
if not 'SOMATICSEQOUT' in globals():
    SOMATICSEQOUT = OUTDIR + 'variants/somatic_seq/'
if not 'STRELKAIN' in globals():
    STRELKAIN = BASERECALIBRATIONOUT
if not 'STRELKAOUT' in globals():
    STRELKAOUT = OUTDIR + 'variants/strelka/'
if not 'STRELKAFILTERIN' in globals():
    STRELKAFILTERIN = STRELKAOUT
if not 'STRELKAFILTEROUT' in globals():
    STRELKAFILTEROUT = OUTDIR + 'variants/strelka/filtered/'

# Definition of some constantly used lists

SAMPLENAMES=getSampleNames()
SINGLEFASTQFILES=getSingleFastqFiles(SAMPLENAMES)
PAIREDFASTQFILES=getPairedFastqFiles(SAMPLENAMES)
PAIREDFASTQFILESWITHOUTR=getPairedFastqFilesWithoutR(SAMPLENAMES)

# Include the rules
include: "../common/clip_trim/clip_trim_snake.py"
include: "../common/align/align_snake.py"
include: "../common/bam_mod/bam_mod_snake.py"
include: "../common/stats/stats_snake.py"
include: "../common/copy_number/copy_number_snake.py"
include: "../common/variants/variants_snake.py"
include: "../common/variants/somatic_seq_snake.py"
include: "../common/variants/variant_mod_snake.py"

# This rule defines which files should be created
rule wes:
    input:
        [file.replace('.fastq.gz', '_fastqc.html') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '_fastqc.html').replace(FASTQDIR, TRIMMOMATICOUT) for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '_fastqc.html').replace(FASTQDIR, TRIMMOMATICOUT).replace('PAIREDEND', 'PAIREDEND/ORPHAN') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '.count') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '.count').replace(FASTQDIR, TRIMMOMATICOUT) for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '.count').replace(FASTQDIR, TRIMMOMATICOUT).replace('PAIREDEND', 'PAIREDEND/ORPHAN') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        expand(MERGEBAMSOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf', sample = SAMPLENAMES),
        #expand(JOINTSNVMIX2_075_OUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        #expand(JOINTSNVMIX2OUT + '{tumorNormalMatching}_classify.jsm', tumorNormalMatching = getNormalTumorFiles()),
        #expand(VARSCANSOMATICOUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        #expand(MUTECT1OUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        #expand(STRELKAOUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        #expand(DEEPSNVOUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        #expand(RANKCOMBINEOUT + '{tumorNormalMatching}.txt', tumorNormalMatching = getNormalTumorFiles()),
        #expand(SOMATICSEQOUT + '{type}.{tumorNormalMatching}.somaticseq.dbSNP.cosmic.snpEff.vcf', type = ['snp', 'indel'], tumorNormalMatching = getNormalTumorFiles()),
        #HAPLOTYPECALLEROUT + 'combined_dist.pdf',
        #expand(FACETSOUT + '{tumorNormalMatching}.cn', tumorNormalMatching = getNormalTumorFiles()),
        expand(VARSCANCNVOUT + '{tumorNormalMatching}.copynumber', tumorNormalMatching = getNormalTumorFiles()),
        expand(GATKVARIANTCOMBINEOUT + '{tumorNormalMatching}.combined.dbSNP.cosmic.snpEff.vcf', tumorNormalMatching = getNormalTumorFiles())
    output:
        OUTDIR + 'complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log',
        mem = '1000',
        scratch = '1000',
        time = '1'
    benchmark:
        OUTDIR + 'complete.txt.benchmark'
    shell:
        'date > {output}'
