# This function returns the TSV file connected to a given FASTQ file
def getTSV(wildcards):
    localFastq =''
    if wildcards.fastq[-3:] == '_R1' or wildcards.fastq[-3:] == '-R1' or wildcards.fastq[-3:] == '_R2' or wildcards.fastq[-3:] == '-R2':
        localFastq = wildcards.fastq[0:-3]
    else:
        localFastq = wildcards.fastq
    return FASTQDIR + wildcards.sample + '/PAIREDEND/' + localFastq + '.tsv'

# This function extracts the information to create the read-group identifier
# given to the read mappers
def createReadGroup(wildcards):
    out = [] # list containing flowcellID, lane, library, platform in that order
    localFastq =''
    if wildcards.fastq[-3:] == '_R1' or wildcards.fastq[-3:] == '-R1' or wildcards.fastq[-3:] == '_R2' or wildcards.fastq[-3:] == '-R2':
        localFastq = wildcards.fastq[0:-3]
    else:
        localFastq = wildcards.fastq
    tsv = open(FASTQDIR+ wildcards.sample + '/PAIREDEND/' + localFastq + '.tsv')
    flowcellID = ''
    lane = ''
    library = ''
    platform = ''
    for line in tsv:
        if line.strip().startswith('RUN_NAME_FOLDER'):
            if flowcellID == '':
                flowcellID = line.strip().split('\t')[1]
        elif line.strip().startswith('LANE_NUMBER'):
            if lane == '':
                lane = line.strip().split('\t')[1]
        elif line.strip().startswith('SAMPLE_CODE'):
            if library == '':
                library = line.strip().split('\t')[1]
        elif line.strip().startswith('SAMPLE_TYPE'):
            if platform == '':
                platform = line.strip().split('\t')[1]
    tsv.close()
    out.append(flowcellID)
    out.append(lane)
    out.append(library)
    out.append(platform)
    return out

# This function creates the read-group entry for bowtie2
def createReadGroupBowtie2(wildcards):
    values = createReadGroup(wildcards)
    return '--rg-id ' + wildcards.sample + '.' + values[0] + '.' + values[1] + ' --rg LB:' + values[2] + ' --rg PL:' + values[3] + ' --rg PU:' + values[0] + '.' + values[1] + '.' + values[2] + ' --rg SM:' + wildcards.sample

# This function creates the read-group entry for yara
def createReadGroupYara(wildcards):
    values = createReadGroup(wildcards)
    return '\'@RG\\tID:' + wildcards.sample + '.' + values[0] + '.' + values[1] + '\\tLB:' + values[2] + '\\tPL:' + values[3] + '\\tPU:' + values[0] + '.' + values[1] + '.' + values[2]  + '\\tSM:' + wildcards.sample + '\''

# This function creates the read-group entry for bwa
def createReadGroupBwa(wildcards):
    values = createReadGroup(wildcards)
    return '\'@RG\\tID:' + wildcards.sample + '.' + values[0] + '.' + values[1] + '\\tLB:' + values[2] + '\\tPL:' + values[3] + '\\tPU:' + values[0] + '.' + values[1] + '.' + values[2]  + '\\tSM:' + wildcards.sample + '\''


# This rule aligns unpaired reads using bowtie2
if not 'BOWTIEIN' in globals():
    BOWTIEIN = CLIPTRIMOUT
if not 'BOWTIEOUT' in globals():
    BOWTIEOUT = OUTDIR + 'bowtie2_out/'
rule alignUnpairedFilesBowtie2:
    input:
        fastq = BOWTIEIN + '{sample}/PAIREDEND/ORPHAN/{fastq}.fastq.gz',
        index1 = config['resources'][ORGANISM]['bowtie2Index'] + '.1.bt2',
        index2 = config['resources'][ORGANISM]['bowtie2Index'] + '.2.bt2',
        index3 = config['resources'][ORGANISM]['bowtie2Index'] + '.3.bt2',
        index4 = config['resources'][ORGANISM]['bowtie2Index'] + '.4.bt2',
        index5 = config['resources'][ORGANISM]['bowtie2Index'] + '.rev.1.bt2',
        index6 = config['resources'][ORGANISM]['bowtie2Index'] + '.rev.2.bt2',
        tsv = getTSV,
    output:
        bam=temp(BOWTIEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam'),
    params:
        sensitivity = config['tools']['bowtie2']['single']['sensitivity'],
        k = config['tools']['bowtie2']['single']['k'],
        phred = config['tools']['bowtie2']['single']['phred'],
        lsfoutfile = BOWTIEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.lsfout.log',
        lsferrfile = BOWTIEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bowtie2']['scratch'],
        mem = config['tools']['bowtie2']['mem'],
        time = config['tools']['bowtie2']['time'],
        index = config['resources'][ORGANISM]['bowtie2Index'],
        rg = createReadGroupBowtie2
    benchmark:
        BOWTIEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.benchmark'
    threads:
        config['tools']['bowtie2']['threads']
    log:
        BOWTIEOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.log'
    shell:
        ('{config[tools][bowtie2][call]} ' +
        '-x {params.index} ' +
        '-U {input.fastq} ' +
        '-k {params.k} ' +
        '{params.phred} ' +
        '-p {threads} ' +
        '{params.sensitivity} ' +
        '{params.rg} ' +
        '2> {log} | ' +
        '{config[tools][samtools][call]} view -bhS - >{output.bam}')

# This rule aligns paired reads using bowtie2
rule alignPairedFilesBowtie2:
    input:
        fastqR1 = BOWTIEIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        fastqR2 = BOWTIEIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        index1 = config['resources'][ORGANISM]['bowtie2Index'] + '.1.bt2',
        index2 = config['resources'][ORGANISM]['bowtie2Index'] + '.2.bt2',
        index3 = config['resources'][ORGANISM]['bowtie2Index'] + '.3.bt2',
        index4 = config['resources'][ORGANISM]['bowtie2Index'] + '.4.bt2',
        index5 = config['resources'][ORGANISM]['bowtie2Index'] + '.rev.1.bt2',
        index6 = config['resources'][ORGANISM]['bowtie2Index'] + '.rev.2.bt2',
        tsv = getTSV,
    output:
        bam=temp(BOWTIEOUT + '{sample}/PAIREDEND/{fastq}.bam'),
    params:
        sensitivity = config['tools']['bowtie2']['single']['sensitivity'],
        k = config['tools']['bowtie2']['single']['k'],
        phred = config['tools']['bowtie2']['single']['phred'],
        lsfoutfile = BOWTIEOUT + '{sample}/PAIREDEND/{fastq}.bam.lsfout.log',
        lsferrfile = BOWTIEOUT + '{sample}/PAIREDEND/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bowtie2']['scratch'],
        mem = config['tools']['bowtie2']['mem'],
        time = config['tools']['bowtie2']['time'],
        index = config['resources'][ORGANISM]['bowtie2Index'],
        rg = createReadGroupBowtie2
    benchmark:
        BOWTIEOUT + '{sample}/PAIREDEND/{fastq}.bam.benchmark'
    threads:
        config['tools']['bowtie2']['threads']
    log:
        BOWTIEOUT + '{sample}/PAIREDEND/{fastq}.bam.log'
    shell:
        ('{config[tools][bowtie2][call]} ' +
        '-x {params.index} ' +
        '-1 {input.fastqR1} ' +
        '-2 {input.fastqR2} ' +
        '-k {params.k} ' +
        '{params.phred} ' +
        '-p {threads} ' +
        '{params.sensitivity} ' +
        '{params.rg} ' +
        '2> {log} | ' +
        '{config[tools][samtools][call]} view -bhS - >{output.bam}')

# This rule aligns unpaired reads using bwa-mem
if not 'BWAIN' in globals():
    BWAIN = CLIPTRIMOUT
if not 'BWAOUT' in globals():
    BWAOUT = OUTDIR + 'bwa/'
rule alignUnpairedFilesBwa:
    input:
        fastq = BWAIN + '{sample}/PAIREDEND/ORPHAN/{fastq}.fastq.gz',
        index1 = config['resources'][ORGANISM]['bwaIndex'] + '.amb',
        index2 = config['resources'][ORGANISM]['bwaIndex'] + '.ann',
        index3 = config['resources'][ORGANISM]['bwaIndex'] + '.bwt',
        index4 = config['resources'][ORGANISM]['bwaIndex'] + '.pac',
        index5 = config['resources'][ORGANISM]['bwaIndex'] + '.sa',
        tsv = getTSV
    output:
        bam=temp(BWAOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam')
    params:
        lsfoutfile = BWAOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.lsfout.log',
        lsferrfile = BWAOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bwa']['mem']['scratch'],
        mem = config['tools']['bwa']['mem']['memory'],
        time = config['tools']['bwa']['mem']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        params = config['tools']['bwa']['mem']['params'],
        rg = createReadGroupBwa
    benchmark:
        BWAOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.benchmark'
    threads:
        config['tools']['bwa']['mem']['threads']
    log:
        BWAOUT + '{sample}/PAIREDEND/ORPHAN/{fastq}.bam.log'
    shell:
        ('{config[tools][bwa][mem][call]} ' +
        '{params.params} ' +
        '-R {params.rg} ' +
        '-t {threads} ' +
        '{params.index} ' +
        '{input.fastq} ' + 
        '2>{log} | {config[tools][samtools][call]} view -bhS - > {output.bam}')

# This rule aligns unpaired reads using bwa-mem
rule alignPairedFilesBwaMem:
    input:
        fastqR1 = BWAIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        fastqR2 = BWAIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        index1 = config['resources'][ORGANISM]['bwaIndex'] + '.amb',
        index2 = config['resources'][ORGANISM]['bwaIndex'] + '.ann',
        index3 = config['resources'][ORGANISM]['bwaIndex'] + '.bwt',
        index4 = config['resources'][ORGANISM]['bwaIndex'] + '.pac',
        index5 = config['resources'][ORGANISM]['bwaIndex'] + '.sa',
        tsv = getTSV
    output:
        bam=temp(BWAOUT + '{sample}/PAIREDEND/{fastq}.bam')
    params:
        lsfoutfile = BWAOUT + '{sample}/PAIREDEND/{fastq}.bam.lsfout.log',
        lsferrfile = BWAOUT + '{sample}/PAIREDEND/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bwa']['mem']['scratch'],
        mem = config['tools']['bwa']['mem']['memory'],
        time = config['tools']['bwa']['mem']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        params = config['tools']['bwa']['mem']['params'],
        rg = createReadGroupBwa
    benchmark:
        BWAOUT + '{sample}/PAIREDEND/{fastq}.bam.benchmark'
    threads:
        config['tools']['bwa']['mem']['threads']
    log:
        BWAOUT + '{sample}/PAIREDEND/{fastq}.bam.log'
    shell:
        ('{config[tools][bwa][mem][call]} ' +
        '{params.params} ' +
        '-R {params.rg} ' +
        '-t {threads} ' +
        '{params.index} ' +
        '{input.fastqR1} ' +
        '{input.fastqR2} ' +
        '2>{log}| {config[tools][samtools][call]} view -bhS - > {output.bam}')

if not 'BWAALNIN' in globals():
    BWAALNIN = CLIPTRIMOUT
if not 'BWAALNOUT' in globals():
    BWAALNOUT = OUTDIR + 'bwa_aln/'
rule alignPairedFilesBwaAln:
    input:
        fastq = BWAALNIN + '{sample}/PAIREDEND/{fastq}.fastq.gz',
        index1 = config['resources'][ORGANISM]['bwaIndex'] + '.amb',
        index2 = config['resources'][ORGANISM]['bwaIndex'] + '.ann',
        index3 = config['resources'][ORGANISM]['bwaIndex'] + '.bwt',
        index4 = config['resources'][ORGANISM]['bwaIndex'] + '.pac',
        index5 = config['resources'][ORGANISM]['bwaIndex'] + '.sa',
    output:
        sai = temp(BWAALNOUT + '{sample}/PAIREDEND/{fastq}.sai')
    params:
        lsfoutfile = BWAALNOUT + '{sample}/PAIREDEND/{fastq}.sai.lsfout.log',
        lsferrfile = BWAALNOUT + '{sample}/PAIREDEND/{fastq}.sai.lsferr.log',
        scratch = config['tools']['bwa']['aln']['scratch'],
        mem = config['tools']['bwa']['aln']['memory'],
        time = config['tools']['bwa']['aln']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        trimQual = config['tools']['bwa']['aln']['trimQual'],
    benchmark:
        BWAALNOUT + '{sample}/PAIREDEND/{fastq}.sai.benchmark'
    threads:
        config['tools']['bwa']['aln']['threads']
    log:
        BWAALNOUT + '{sample}/PAIREDEND/{fastq}.sai.log'
    shell:
        ('{config[tools][bwa][aln][call]} ' +
        '{params.trimQual} ' +
        '-t {threads} ' +
        '{params.index} ' +
        '{input.fastq} ' +
        '2>{log}| > {output.sai}')

rule alignPairedFilesBwaSampe:
    input:
        sai1 = BWAALNOUT + '{sample}/PAIREDEND/{fastq}_R1.sai',
        fastq_R1 = BWAALNIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        sai2 = BWAALNOUT + '{sample}/PAIREDEND/{fastq}_R2.sai',
        fastq_R2 = BWAALNIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
        tsv = getTSV
    output:
        bam=temp(BWAALNOUT + '{sample}/PAIREDEND/{fastq}.bam')
    params:
        lsfoutfile = BWAALNOUT + '{sample}/PAIREDEND/{fastq}.bam.lsfout.log',
        lsferrfile = BWAALNOUT + '{sample}/PAIREDEND/{fastq}.bam.lsferr.log',
        scratch = config['tools']['bwa']['sampe']['scratch'],
        mem = config['tools']['bwa']['sampe']['mem'],
        time = config['tools']['bwa']['sampe']['time'],
        index = config['resources'][ORGANISM]['bwaIndex'],
        rg = createReadGroupBwa
    benchmark:
        BWAALNOUT + '{sample}/PAIREDEND/{fastq}.bam.benchmark'
    threads:
        config['tools']['bwa']['sampe']['threads']
    log:
        BWAALNOUT + '{sample}/PAIREDEND/{fastq}.bam.log'
    shell:
        ('{config[tools][bwa][sampe][call]} ' +
        '-r {params.rg} ' +
        '{params.index} ' +
        '{input.sai1} ' +
        '{input.sai2} ' +
        '{input.fastq_R1} ' +
        '{input.fastq_R2} ' +
        '2>{log}| {config[tools][samtools][call]} view -bhS - > {output.bam}')

# This rule aligns unpaired reads using yara
if not 'YARAIN' in globals():
    YARAIN = OUTDIR + '.yara_in'
if not 'YARAOUT' in globals():
    YARAOUT = OUTDIR + '.yara_out'
rule alignUnpairedFilesYara:
    input:
        fastq = YARAIN + '{sample}/PAIREDEND/{fastq}_UNPAIRED.fastq.gz',
        index1 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drp',
        index2 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drs',
        index3 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drv',
        index4 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.pst',
        index5 = config['resources'][ORGANISM]['yaraIndex'] + '.rid.concat',
        index6 = config['resources'][ORGANISM]['yaraIndex'] + '.rid.limits',
        index7 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.ind',
        index8 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.len',
        index9 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.val',
        index10 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.concat',
        index11 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.limits',
        index12 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.size',
        tsv = getTSV
    output:
        bam=temp(YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam'),
    params:
        lsfoutfile = YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.lsfout.log',
        lsferrfile = YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.lsferr.log',
        scratch = config['tools']['yara']['scratch'],
        mem = config['tools']['yara']['mem'],
        time = config['tools']['yara']['time'],
        index = config['resources'][ORGANISM]['yaraIndex'],
        rg = createReadGroupYara
    benchmark:
        YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.benchmark'
    threads:
        config['tools']['yara']['threads']
    log:
        YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.log'
    shell:
        ('{config[tools][yara][call]} ' +
        '-rg {params.rg} ' +
        '-o {output.bam} ' +
        '-e config[tools][yara][paired][error-rate] ' +
        '-s config[tools][yara][paired][strata-rate] ' +
        '-t {threads} ' +
        '{input.fastq} ' +
        '2> {log}')

# This rule aligns paired reads using yara
rule alignPairedFilesYara:
    input:
        fastqR1 = YARAIN + '{sample}/PAIREDEND/{fastq}_R1_PAIRED.fastq.gz',
        fastqR2 = YARAIN + '{sample}/PAIREDEND/{fastq}_R2_PAIRED.fastq.gz',
        index1 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drp',
        index2 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drs',
        index3 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.drv',
        index4 = config['resources'][ORGANISM]['yaraIndex'] + '.lf.pst',
        index5 = config['resources'][ORGANISM]['yaraIndex'] + '.rid.concat',
        index6 = config['resources'][ORGANISM]['yaraIndex'] + '.rid.limits',
        index7 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.ind',
        index8 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.len',
        index9 = config['resources'][ORGANISM]['yaraIndex'] + '.sa.val',
        index10 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.concat',
        index11 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.limits',
        index12 = config['resources'][ORGANISM]['yaraIndex'] + '.txt.size',
        tsv = getTSV
    output:
        bam=temp(YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam'),
    params:
        lsfoutfile = YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.lsfout.log',
        lsferrfile = YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.lsferr.log',
        scratch = config['tools']['yara']['scratch'],
        mem = config['tools']['yara']['mem'],
        time = config['tools']['yara']['time'],
        index = config['resources'][ORGANISM]['yaraIndex'],
        rg = createReadGroupYara
    benchmark:
        YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.benchmark'
    threads:
        config['tools']['yara']['threads']
    log:
        YARAOUT + '{sample}/PAIREDEND/{fastq}_PAIRED.bam.log'
    shell:
        ('{config[tools][yara][call]} ' +
        '-rg {params.rg} ' +
        '-o {output.bam} ' +
        '-e config[tools][yara][paired][error-rate] ' +
        '-s config[tools][yara][paired][strata-rate] ' +
        '-t {threads} ' +
        '{input.fastqR1} ' +
        '{input.fastqR2} ' +
        '2> {log}')

# This rule aligns unpaired reads using bwa-mem
if not 'SOAPIN' in globals():
    SOAPIN = OUTDIR + '.yara_in'
if not 'SOAPOUT' in globals():
    SOAPOUT = OUTDIR + '.soap_out'
rule alignPairedFilesSoap:
    input:
        fastqR1 = SOAPIN + '{sample}/PAIREDEND/{fastq}_R1.fastq',
        fastqR2 = SOAPIN + '{sample}/PAIREDEND/{fastq}_R2.fastq',
        index = config['resources'][ORGANISM]['soapIndex'] + '.bwt'
    output:
        txt=temp(SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt'),
        txtSE=temp(SOAPOUT + '{sample}/PAIREDEND/{fastq}_SE.txt')
    params:
        lsfoutfile = SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt.lsfout.log',
        lsferrfile = SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt.lsferr.log',
        scratch = config['tools']['soap']['scratch'],
        mem = config['tools']['soap']['memory'],
        time = config['tools']['soap']['time'],
        params = config['tools']['soap']['params'],
        index = config['resources'][ORGANISM]['soapIndex'],
        rg = createReadGroupBwa
    benchmark:
        SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt.benchmark'
    threads:
        config['tools']['soap']['threads']
    log:
        SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt.log'
    shell:
        ('{config[tools][soap][call]} ' +
        '-a {input.fastqR1} ' +
        '-b {input.fastqR2} ' +
        '-D {params.index} ' +
        '-o {output.txt} ' + 
        '-2 {output.txtSE} ' +
        '-p {threads} ' +
        '{params.params} ')

# This rule aligns unpaired reads using bwa-mem
rule soap2sam:
    input:
        txt = SOAPOUT + '{sample}/PAIREDEND/{fastq}.txt',
        fai = config['resources'][ORGANISM]['referenceFai']
    output:
        bam=temp(SOAPOUT + '{sample}/PAIREDEND/{fastq}.bam')
    params:
        lsfoutfile = SOAPOUT + '{sample}/PAIREDEND/{fastq}.bam.lsfout.log',
        lsferrfile = SOAPOUT + '{sample}/PAIREDEND/{fastq}.bam.lsferr.log',
        scratch = config['tools']['soap2sam']['scratch'],
        mem = config['tools']['soap2sam']['memory'],
        time = config['tools']['soap2sam']['time'],
        rg = createReadGroupBwa
    benchmark:
        SOAPOUT + '{sample}/PAIREDEND/{fastq}.bam.benchmark'
    threads:
        config['tools']['soap2sam']['threads']
    log:
        SOAPOUT + '{sample}/PAIREDEND/{fastq}.bam.log'
    shell:
        ('{config[tools][soap2sam][call]} ' +
        '-p ' +
        '{input.txt} '
        '| {config[tools][samtools][call]} view -bt {input.fai} - >{output.bam}')

