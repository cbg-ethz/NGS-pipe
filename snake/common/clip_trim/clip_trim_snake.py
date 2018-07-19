# This rule clips bases with too low quality from the reads
# and removes adapter contamination in the reads
if not 'TRIMMOMATICIN' in globals():
    TRIMMOMATICIN = FASTQDIR
if not 'TRIMMOMATICOUT' in globals():
    TRIMMOMATICOUT = OUTDIR + 'cliptrim/'
rule trimmomatic_paired:
    input:
        forward = TRIMMOMATICIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        reverse = TRIMMOMATICIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        adapter = config['resources']['general']['sequencingAdapter']
    output:
        forwardP = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz'),
        forwardUP = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}_R1.fastq.gz'),
        reverseP = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'),
        reverseUP = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}_R2.fastq.gz'),
        trimlog = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log.gz'
    params:
        trimlog = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log'),
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
        lsfoutfile = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.lsfout.log',
        lsferrfile = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.lsferr.log',
        scratch = config['tools']['trimmomatic']['scratch'],
        mem = config['tools']['trimmomatic']['mem'],
        time = config['tools']['trimmomatic']['time']
    benchmark:
        TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.benchmark'
    threads:
        config['tools']['trimmomatic']['paired']['threads']
    log:
        stdoutlog = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.stdout.log'
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

rule trimmomatic_single:
    input:
        fastq = TRIMMOMATICIN + '{sample}/SINGLEEND/{fastq}.fastq.gz',
        adapter = config['trimmomatic']['single']['rnaadapterfile']
    output:
        temp(TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}.fastq.gz')
    params:
        slidingwindow = config['trimmomatic']['single']['slidingwindow'],
        phred = config['trimmomatic']['single']['phred'],
        mode = config['trimmomatic']['single']['mode'],
        minQual = config['trimmomatic']['single']['minQual'],
        seedmismatches = config['trimmomatic']['single']['seedmismatches'],
        palindrom = config['trimmomatic']['single']['palindrom'],
        score = config['trimmomatic']['single']['score'],
        minlen = config['trimmomatic']['single']['minlen'],
        lsfoutfile = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsfout.log',
        lsferrfile = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsferr.log',
        scratch = config['trimmomatic']['scratch'],
        mem = config['trimmomatic']['mem'],
        time = config['trimmomatic']['time']
    benchmark:
        TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}.benchmark'
    threads:
        config['tools']['trimmomatic']['single']['threads']
    log:
        trimlog = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.log',
        stdoutlog = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.stdout.log'
    shell:
        ('{config[trimmomatic][call]} ' + 
        '{params.mode} ' + 
        '{params.phred} ' + 
        '-threads {threads} ' + 
        '{input.fastq} ' + 
        '{output} ' + 
        'ILLUMINACLIP:{input.adapter}:{params.seedmismatches}:{params.palindrom}:{params.score} ' + 
        'SLIDINGWINDOW:{params.slidingwindow}:{params.minQual} ' + 
        'LEADING:{params.minQual} ' + 
        'TRAILING:{params.minQual} ' + 
        'MINLEN:{params.minlen} ' + 
        '2> {log.stdoutlog}')

if not 'SEQPURGEIN' in globals():
    SEQPURGEIN = FASTQDIR
if not 'SEQPURGEOUT' in globals():
    SEQPURGEOUT = OUTDIR + 'seqpurge/'
rule seqpurge_paired:
    input:
        in1 = SEQPURGEIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        in2 = SEQPURGEIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'
    output:
        out1 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz'),
        out2 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz'),
        out3 = temp(SEQPURGEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.fastq.gz'),
    params:
        lsfoutfile = SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.lsfout.log',
        lsferrfile = SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.lsferr.log',
        scratch = config['tools']['seqpurge']['scratch'],
        mem = config['tools']['seqpurge']['mem'],
        time = config['tools']['seqpurge']['time'],
        a1 = config['tools']['seqpurge']['a1'],
        a2 = config['tools']['seqpurge']['a2'],
        params = config['tools']['seqpurge']['params']
    benchmark:
        SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.benchmark'
    threads:
        config['tools']['seqpurge']['threads']
    log:
        SEQPURGEOUT + '{sample}/PAIREDEND/{fastq}.log'
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
