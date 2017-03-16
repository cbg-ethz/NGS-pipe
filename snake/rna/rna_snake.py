

OUTFOLDER = ROOTFOLDER + 'out/'
CLIPTRIMOUT = OUTFOLDER + 'cliptrim1/'
ALIGNOUT = OUTFOLDER + 'align2/'
TMP = OUTFOLDER + 'tmp/'
STATSOUT = OUTFOLDER + 'stats4/'
FEATURECOUNTSOUT = OUTFOLDER + 'featurecounts3/'
'''
Collecting file names
'''
PAIREDFASTQFILES = [file.replace(FASTQFOLDER, '').replace('.fastq.gz','')for file in glob.glob(FASTQFOLDER + '*/PAIREDEND/*.fastq.gz')]
SINGLEFASTQFILES = [file.replace(FASTQFOLDER, '').replace('.fastq.gz','')for file in glob.glob(FASTQFOLDER + '*/SINGLEEND/*.fastq.gz')]
SAMPLENAMES = [samplename.replace(FASTQFOLDER,'').replace('/','')for samplename in glob.glob(FASTQFOLDER + '*/')]
PAIREDENDSAMPLENAMES = [token for token in set([samplename.split('/')[0].strip() for samplename in PAIREDFASTQFILES])]
SINGLEENDSAMPLENAMES = [token for token in set([samplename.split('/')[0].strip() for samplename in SINGLEFASTQFILES])]
MAPPINGFILES = glob.glob(FASTQFOLDER + '/*/mapping.tsv')




'''
methods
'''

def getFastqsAfterTrimming(wildcards):
    '''
        get the fastqs after trimming
    '''
    paired = expand(CLIPTRIMOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    single = expand(CLIPTRIMOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    all = []
    for i in single:
        all.append(i)
    for i in paired:
        all.append(i)
    return all


def getFastQForSampleSingleEnd(wildcards):
    sample = wildcards.sample
    singles = expand(CLIPTRIMOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    wc2single = []
    for single in singles:
        if '/' + sample+ '/' in single:
            wc2single.append(single)
    return wc2single
def getFastQForSampleSingleEndStarString(wildcards):
    fastqs = getFastQForSampleSingleEnd(wildcards)
    s = ''
    for fastq in fastqs:
        s = s + ',' + fastq
    return s[1:]

def getFastQForSamplePairedEnd(wildcards):
    sample = wildcards.sample
    paireds = expand(CLIPTRIMOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    wc2paired = []
    for paired in paireds:
        if '/' + sample+ '/' in paired:
            wc2paired.append(paired)
    return wc2paired

def getFastQForSamplePairedEndStarString(wildcards):
    fastqs = getFastQForSamplePairedEnd(wildcards)
    left = []
    right = []
    for fastq in fastqs:
        if '_R1_PAIRED' in fastq:
            left.append(fastq)
            continue
        if '_R2_PAIRED' in fastq:
            right.append(fastq)
            continue
        raise ValueError('There re fastqs that do not fit the scheme for this sample: ' + fastq)
    if len(left) != len(right):
        raise ValueError('Not all paired end files are matched: ' + str(fastqs))
    left.sort()
    right.sort()

    return ','.join(left) + ' ' + ','.join(right)



def getAllFastqFiles(wildcards):
    infiles = expand(FASTQFOLDER + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    single = expand(CLIPTRIMOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    inpairedfiles = expand(FASTQFOLDER + '{files}.fastq.gz', files=PAIREDFASTQFILES)
    cliptrimpaired = expand(CLIPTRIMOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    all = single + infiles + inpairedfiles + cliptrimpaired
    return all
def getAllFastqCountFiles(wildcards):
    fastqs = getAllFastqFiles(wildcards)
    counts = []
    for fastq in fastqs:
        counts.append(fastq.replace('fastq.gz', 'count'))
    return counts
def getAllFastqcFiles(wildcards):
    fastqs = getAllFastqFiles(wildcards)
    fastqcs = []
    for fastq in fastqs:
        fastqcs.append(fastq.replace('fastq.gz', 'fastqc.done'))
    return fastqcs
def getAllBamFiles(wildcards):
    alignpaired = expand(ALIGNOUT + 'PAIREDEND/{files}.bam', files=PAIREDENDSAMPLENAMES)
    alignsingle = expand(ALIGNOUT + 'SINGLEEND/{files}.bam', files=SINGLEENDSAMPLENAMES)
    align = alignpaired + alignsingle
    return align
def getBamHeaderFiles(wildcards):
    bams = getAllBamFiles(wildcards)
    headers = []
    for bam in bams:
        headers.append(bam.replace('.bam', '.header'))
    return headers
def getFlagstatFiles(wildcards):
    bams = getAllBamFiles(wildcards)
    flagstats = []
    for bam in bams:
        flagstats.append(bam.replace('.bam', '.flagstat'))
    return flagstats
def getAllStatsFilesForASample(wildcards):
    flagstats = getFlagstatFiles(None)
    counts = getAllFastqCountFiles(None)
    sample = wildcards.sample
    out = []
    for flagstat in flagstats:
        if sample in flagstat:
            out.append(flagstat)
    for count in counts:
        if sample in count:
            out.append(count)
    return out
def getAllAlignFiles(wildcards):
    return expand(ALIGNOUT + 'PAIREDEND/{files}.bam', files=PAIREDENDSAMPLENAMES) + expand(ALIGNOUT + 'SINGLEEND/{files}.bam', files=SINGLEENDSAMPLENAMES)
def getAllFcntFiles(wildcards):
    return expand(FEATURECOUNTSOUT + 'PAIREDEND/{sample}.fcnt', sample=PAIREDENDSAMPLENAMES) + expand(FEATURECOUNTSOUT + 'SINGLEEND/{sample}.fcnt', sample=SINGLEENDSAMPLENAMES)
def getAllFcnt2HtseqFiles(wildcards):
    return expand(FEATURECOUNTSOUT + 'PAIREDEND/{sample}.fcnt2htseq', sample=PAIREDENDSAMPLENAMES) + expand(FEATURECOUNTSOUT + 'SINGLEEND/{sample}.fcnt2htseq', sample=SINGLEENDSAMPLENAMES)

'''
rules

'''
rule all:
    input:
        getBamHeaderFiles,
        getAllFastqcFiles,
        getFlagstatFiles,
        getAllFastqCountFiles,
        getFastqsAfterTrimming,
        getAllAlignFiles,
        getAllFcntFiles,
        getAllFcnt2HtseqFiles
    output:
        touch('complete.txt')


rule trimSingle:
    '''
        trim the single end reads with trimmomatic
        '''
    input:
        fastq = FASTQFOLDER + '{sample}/SINGLEEND/{fastq}.fastq.gz',
        adapter = config['trimmomatic']['single']['rnaadapterfile']
    output:
        temp(CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}.fastq.gz')
    params:
        slidingwindow = config['trimmomatic']['single']['slidingwindow'],
        phred = config['trimmomatic']['single']['phred'],
        mode = config['trimmomatic']['single']['mode'],
        minQual = config['trimmomatic']['single']['minQual'],
        seedmismatches = config['trimmomatic']['single']['seedmismatches'],
        palindrom = config['trimmomatic']['single']['palindrom'],
        score = config['trimmomatic']['single']['score'],
        minlen = config['trimmomatic']['single']['minlen'],
        lsfoutfile = CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsfout.log',
        lsferrfile = CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsferr.log',
        scratch = config['trimmomatic']['scratch'],
        mem = config['trimmomatic']['mem'],
        time = config['trimmomatic']['time']
    benchmark:
        CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}.benchmark'
    threads:
        int(config['trimmomatic']['single']['threads'])
    log:
        trimlog = CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.log',
        stdoutlog = CLIPTRIMOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.stdout.log'
    shell:
        '{config[trimmomatic][call]} {params.mode} {params.phred} -threads {threads} {input.fastq} {output} ILLUMINACLIP:{input.adapter}:{params.seedmismatches}:{params.palindrom}:{params.score} SLIDINGWINDOW:{params.slidingwindow}:{params.minQual} LEADING:{params.minQual} TRAILING:{params.minQual} MINLEN:{params.minlen} 2> {log.stdoutlog}'

rule clipTrimPaired:
    input:
        forward = FASTQFOLDER + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        reverse = FASTQFOLDER + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        adapter = config['trimmomatic']['paired']['rnaadapterfile']
    output:
        forwardP = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R1_PAIRED.fastq.gz',
        forwardUP = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R1_UNPAIRED.fastq.gz',
        reverseP = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R2_PAIRED.fastq.gz',
        reverseUP = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R2_UNPAIRED.fastq.gz',
        trimlog = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log'
    params:
        slidingwindow = config['trimmomatic']['paired']['slidingwindow'],
        phred = config['trimmomatic']['paired']['phred'],
        mode = config['trimmomatic']['paired']['mode'],
        minQual = config['trimmomatic']['paired']['minQual'],
        seedmismatches = config['trimmomatic']['paired']['seedmismatches'],
        palindrom = config['trimmomatic']['paired']['palindrom'],
        score = config['trimmomatic']['paired']['score'],
        minlen = config['trimmomatic']['paired']['minlen'],
        min_adapt_len = config['trimmomatic']['paired']['min_adapt_len'],
        keep_both = config['trimmomatic']['paired']['keep_both'],
        lsfoutfile = CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}_clipTrim.lsfout.log',
        lsferrfile = CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}_clipTrim.lsferr.log',
        scratch = config['trimmomatic']['scratch'],
        mem = config['trimmomatic']['mem'],
        time = config['trimmomatic']['time']
    benchmark:
        CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}.benchmark'
    threads:
        int(config['trimmomatic']['paired']['threads'])
    log:
        stdoutlog = CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.stdout.log'
    shell:
        ('{config[trimmomatic][call]} ' +
        '{params.mode} ' +
        '{params.phred} ' +
        '-threads {threads} ' +
        '-trimlog {output.trimlog} ' +
        '{input.forward} ' +
        '{input.reverse} ' +
        '{output.forwardP} ' +
        '{output.forwardUP} ' +
        '{output.reverseP} ' +
        '{output.reverseUP} ' +
        'ILLUMINACLIP:{input.adapter}:{params.seedmismatches}:{params.palindrom}:{params.score}:{params.min_adapt_len}:{params.keep_both} ' +
        'SLIDINGWINDOW:{params.slidingwindow}:{params.minQual} ' +
        'LEADING:{params.minQual} ' +
        'TRAILING:{params.minQual} ' +
        'MINLEN:{params.minlen} ' +
        '2> {log.stdoutlog}')

rule alignSingleFileStar:
    '''
        align all fastq files that are single
        '''
    input:
        fastq = getFastQForSampleSingleEnd,
        index = STARINDEXFOLDER
    output:
        ALIGNOUT + 'SINGLEEND/{sample}.bam'
    params:
        fastq = getFastQForSampleSingleEndStarString,
        readFilesCommand = config['star']['readFilesCommand'],
        outFilterMultimapNmax = config['star']['outFilterMultimapNmax'],
        outStd = config['star']['outStd'],
        outSamUnmapped = config['star']['outSamUnmapped'],
        outFilterMismatchNMAX = config['star']['outFilterMismatchNMAX'],
        lsfoutfile = ALIGNOUT + 'SINGLEEND/{sample}_star.lsfout.log',
        lsferrfile = ALIGNOUT + 'SINGLEEND/{sample}_star.lsferr.log',
        scratch = config['star']['scratch'],
        mem = config['star']['mem'],
        time = config['star']['time']
    benchmark:
        ALIGNOUT + 'SINGLEEND/{sample}_star.benchmark'
    threads:
        int(config['star']['threads'])
    log:
        ALIGNOUT + 'SINGLEEND/{sample}_star.log'
    shell:
        '{config[star][call]} --genomeDir {input.index} --readFilesIn {params.fastq} --readFilesCommand {params.readFilesCommand} --runThreadN {threads} --outFileNamePrefix {output} --outFilterMultimapNmax {params.outFilterMultimapNmax} --outStd {params.outStd} --outSAMunmapped {params.outSamUnmapped} --outFilterMismatchNmax {params.outFilterMismatchNMAX} | {config[samtools][call]} view -Shb - > {output}'



rule alignPairedFileStar:
    '''
        align all fastq files that are paired
        '''
    input:
        fastq = getFastQForSamplePairedEnd,
        index = STARINDEXFOLDER
    output:
        ALIGNOUT + 'PAIREDEND/{sample}.bam'
    params:
        fastq = getFastQForSamplePairedEndStarString,
        readFilesCommand = config['star']['readFilesCommand'],
        outFilterMultimapNmax = config['star']['outFilterMultimapNmax'],
        outStd = config['star']['outStd'],
        outSamUnmapped = config['star']['outSamUnmapped'],
        outFilterMismatchNMAX = config['star']['outFilterMismatchNMAX'],
        lsfoutfile = ALIGNOUT + 'PAIREDEND/{sample}_star.lsfout.log',
        lsferrfile = ALIGNOUT + 'PAIREDEND/{sample}_star.lsferr.log',
        scratch = config['star']['scratch'],
        mem = config['star']['mem'],
        time = config['star']['time']
    benchmark:
        ALIGNOUT + 'PAIREDEND/{sample}_star.benchmark'
    threads:
        4
    log:
        ALIGNOUT + 'PAIREDEND/{sample}_star.log'
    shell:
        '{config[star][call]} --genomeDir {input.index} --readFilesIn {params.fastq} --readFilesCommand {params.readFilesCommand} --runThreadN {threads} --outFileNamePrefix {output} --outFilterMultimapNmax {params.outFilterMultimapNmax} --outStd {params.outStd} --outSAMunmapped {params.outSamUnmapped} --outFilterMismatchNmax {params.outFilterMismatchNMAX} | {config[samtools][call]} view -Shb - > {output}'




rule featurecountsSingle:
    input:
        bam = ALIGNOUT + 'SINGLEEND/{sample}.bam',
        gtf = GTFFILE
    output:
        FEATURECOUNTSOUT + 'SINGLEEND/{sample}.fcnt'
    params:
        lsfoutfile = FEATURECOUNTSOUT + 'SINGLEEND/{sample}.lsfout.log',
        lsferrfile = FEATURECOUNTSOUT + 'SINGLEEND/{sample}.lsferr.log',
        scratch = config['featureCounts']['scratch'],
        mem = config['featureCounts']['mem'],
        time = config['featureCounts']['time'],
        parameters = config['featureCounts']['parametersSingle']
    benchmark:
        FEATURECOUNTSOUT + 'SINGLEEND/{sample}.benchmark'
    threads:
        int(config['featureCounts']['threads'])
    shell:
        '{config[featureCounts][call]} -T {threads} {params.parameters} -a {input.gtf} -o {output} {input.bam}'


rule featurecountsPaired:
    input:
        bam = ALIGNOUT + 'PAIREDEND/{sample}.bam',
        gtf = GTFFILE
    output:
        FEATURECOUNTSOUT + 'PAIREDEND/{sample}.fcnt'
    params:
        lsfoutfile = FEATURECOUNTSOUT + 'PAIREDEND/{sample}.lsfout.log',
        lsferrfile = FEATURECOUNTSOUT + 'PAIREDEND/{sample}.lsferr.log',
        scratch = config['featureCounts']['scratch'],
        mem = config['featureCounts']['mem'],
        time = config['featureCounts']['time'],
        parameters = config['featureCounts']['parametersPaired']
    benchmark:
        FEATURECOUNTSOUT + 'PAIREDEND/{sample}.benchmark'
    threads:
        int(config['featureCounts']['threads'])
    shell:
        '{config[featureCounts][call]} -T {threads} {params.parameters} -a {input.gtf} -o {output} {input.bam}'

rule featurecounts2HtSeq:
    input:
        fcnt = FEATURECOUNTSOUT + '{endtype}/{sample}.fcnt'
    output:
        htseq = FEATURECOUNTSOUT + '{endtype}/{sample}.fcnt2htseq'
    params:
        scratch = config['featurecounts2HtSeq']['scratch'],
        mem = config['featurecounts2HtSeq']['mem'],
        time = config['featurecounts2HtSeq']['time'],
        lsfoutfile = FEATURECOUNTSOUT + '{endtype}/{sample}.fcnt2htseq.lsfout.log',
        lsferrfile = FEATURECOUNTSOUT + '{endtype}/{sample}.fcnt2htseq.lsferr.log'
    threads:
        int(config['featurecounts2HtSeq']['threads'])
    shell:
        '{config[featurecounts2HtSeq][call]}  {input.fcnt} {output.htseq}'

rule countReadsInFastq:
    input:
        fastq = '{sample}.fastq.gz'
    output:
        count = '{sample}.count'
    params:
        scratch = config['countReadsInFastq']['scratch'],
        mem = config['countReadsInFastq']['mem'],
        time = config['countReadsInFastq']['time'],
        lsfoutfile = '{sample}.count.lsfout.log',
        lsferrfile = '{sample}.count.lsferr.log'
    threads:
        int(config['countReadsInFastq']['threads'])
    priority:
        2
    shell:
        'zcat {input.fastq} | wc -l | awk \'{{print $1/4}}\' > {output.count}'

rule runFastqcOnFastqs:
    input:
        fastq = '{sample}.fastq.gz'
    output:
        fastqc = touch('{sample}.fastqc.done')
    params:
        scratch = config['fastqc']['scratch'],
        mem = config['fastqc']['mem'],
        time = config['fastqc']['time'],
        lsfoutfile = '{sample}.fastqc.lsfout.log',
        lsferrfile = '{sample}.fastqc.lsferr.log'
    threads:
        int(config['fastqc']['threads'])
    priority:
        2
    shell:
        '{config[fastqc][call]} {input.fastq}'

rule runFlagstat:
    input:
        bam = '{sample}.bam'
    output:
        flagstat = '{sample}.flagstat'
    params:
        mem = config['samtools']['flagstat']['mem'],
        scratch = config['samtools']['flagstat']['scratch'],
        time = config['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.flagstat.lsfout.log',
        lsferrfile = '{sample}.flagstat.lsferr.log'
    threads:
        int(config['samtools']['flagstat']['threads'])
    priority:
        2
    shell:
        '{config[samtools][call]} flagstat {input.bam} > {output.flagstat}'




rule runSamtoolsHeader:
    input:
        bam = '{sample}.bam'
    output:
        header = '{sample}.header'
    params:
        mem = config['samtools']['flagstat']['mem'],
        scratch = config['samtools']['flagstat']['scratch'],
        time = config['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.header.lsfout.log',
        lsferrfile = '{sample}.header.lsferr.log'
    threads:
        int(config['samtools']['flagstat']['threads'])
    priority:
        2
    shell:
        '{config[samtools][call]} view -H {input.bam} > {output.header}'
