# This rule clips bases with too low quality from the reads
# and removes adapter contamination in the reads
if not 'CLIPTRIMIN' in globals():
    CLIPTRIMIN = FASTQDIR
if not 'CLIPTRIMOUT' in globals():
    CLIPTRIMOUT = OUTDIR + 'cliptrim/'
rule clipTrimPaired:
    input:
        forward = CLIPTRIMIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        reverse = CLIPTRIMIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        adapter = config['resources']['general']['sequencingAdapter']
    output:
        forwardP = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz'),
        forwardUP = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}_R1.fastq.gz'),
        reverseP = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'),
        reverseUP = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}_R2.fastq.gz'),
        trimlog = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log.gz')
    params:
        trimlog = temp(CLIPTRIMOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log'),
        slidingwindow = config['tools']['trimmomatic']['paired']['slidingwindow'],
        phred = config['tools']['trimmomatic']['paired']['phred'],
        mode = config['tools']['trimmomatic']['paired']['mode'],
        minQual = config['tools']['trimmomatic']['paired']['minQual'],
        seedmismatches = config['tools']['trimmomatic']['paired']['seedmismatches'],
        palindrom = config['tools']['trimmomatic']['paired']['palindrom'],
        score = config['tools']['trimmomatic']['paired']['score'],
        minlen = config['tools']['trimmomatic']['paired']['minlen'],
        min_adapt_len = config['tools']['trimmomatic']['paired']['min_adapt_len'],
        keep_both = config['tools']['trimmomatic']['paired']['keep_both'],
        lsfoutfile = CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}.fastq.lsfout.log',
        lsferrfile = CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}.fastq.lsferr.log',
        scratch = config['tools']['trimmomatic']['scratch'],
        mem = config['tools']['trimmomatic']['mem'],
        time = config['tools']['trimmomatic']['time']
    benchmark:
        CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}.fastq.gz.benchmark'
    threads:
        config['tools']['trimmomatic']['paired']['threads']
    log:
        stdoutlog = CLIPTRIMOUT + '/{sample}/PAIREDEND/{fastq}.fastq.gz.stdout.log'
    shell:
        ('{config[tools][trimmomatic][call]} ' +
        '{params.mode} ' +
        '{params.phred} ' +
        '-threads {threads} ' +
        '-trimlog {params.trimlog} ' +
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
        '2> {log.stdoutlog} && ' + 
        'gzip {params.trimlog}')

if not 'SEQPURGEIN' in globals():
    SEQPURGEIN = FASTQDIR
if not 'SEQPURGEOUT' in globals():
    SEQPURGEOUT = OUTDIR + 'seqpurge/'
rule SeqPurge:
    input:
        in1 = SEQPURGEIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        in2 = SEQPURGEIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'
    output:
        out1 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz'),
        out2 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'),
        out3 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.fastq.gz'),
    params:
        lsfoutfile = SEQPURGEOUT + '/{sample}/PAIREDEND/{fastq}.fastq.gz.lsfout.log',
        lsferrfile = SEQPURGEOUT + '/{sample}/PAIREDEND/{fastq}.fastq.gz.lsferr.log',
        scratch = config['tools']['seqpurge']['scratch'],
        mem = config['tools']['seqpurge']['mem'],
        time = config['tools']['seqpurge']['time'],
        a1 = config['tools']['seqpurge']['a1'],
        a2 = config['tools']['seqpurge']['a2'],
        params = config['tools']['seqpurge']['params']
    benchmark:
        SEQPURGEOUT + '/{sample}/PAIREDEND/{fastq}.fastq.gz.benchmark'
    threads:
        config['tools']['seqpurge']['threads']
    log:
        SEQPURGEOUT + '/{sample}/PAIREDEND/{fastq}.log'
    shell:
        ('{config[tools][seqpurge][call]} ' +
        '-in1 {input.in1} ' +
        '-in2 {input.in2} ' +
        '-out1 {output.out1} ' +
        '-out2 {output.out2} ' +
        '-out3 {output.out3} ' +
        '-a1 {params.a1} ' +
        '-a2 {params.a2} ' +
        '-summary {log}; ' +
        'touch {output.out3}')
