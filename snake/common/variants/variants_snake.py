
if not 'VARSCANSNPIN' in globals():
    VARSCANSNPIN = MPILEUPOUT
if not 'VARSCANSNPOUT' in globals():
    VARSCANSNPOUT = OUTDIR + 'variants/varscansnp/raw/'
rule varscanSnp:
    input:
        mpileup = VARSCANSNPIN + '{sample}.mpileup'
    output:
        vcf = VARSCANSNPOUT + '{sample}.vcf'
    priority: 25
    params:
        lsfoutfile = VARSCANSNPOUT + '{sample}.vcf.lsfout.log',
        lsferrfile = VARSCANSNPOUT + '{sample}.vcf.lsferr.log'
    threads:
        config['tools']['varscan']['pileup2snp']['threads']
    benchmark:
        VARSCANSNPOUT + '{sample}.vcf.benchmark'
    log:
        VARSCANSNPOUT + '{sample}.vcf.log'
    shell:
        ('{config[tools][varscan][call]} mpileup2snp ' +
        '{input.mpileup} ' +
        '--min-coverage {config[tools][varscan][pileup2snp][minCov]} ' +
        '--min-reads2 {config[tools][varscan][pileup2snp][minReads2]} ' +
        '--min-avg-qual {config[tools][varscan][pileup2snp][minAvgQual]} ' +
        '--min-var-freq {config[tools][varscan][pileup2snp][minVarFreq]} ' +
        '--output-vcf 1 ' +
        '--p-value {config[tools][varscan][pileup2snp][pValue]} ' +
        '2>{log} > {output.vcf} ' +
        '&& touch {output.suc}')

if not 'BCFTOOLSIN' in globals():
    BCFTOOLSIN = MPILEUPOUT
if not 'BCFTOOLSOUT' in globals():
    BCFTOOLSOUT = OUTDIR + 'bcftools/raw/'
rule bcftools:
    input:
        bcf = BCFTOOLSIN + '{sample}.bcf'
    output:
        vcf = BCFTOOLSOUT + '{sample}.vcf'
    priority: 25
    params:
        lsfoutfile = BCFTOOLSOUT + '{sample}.vcf.lsfout.log',
        lsferrfile = BCFTOOLSOUT + '{sample}.vcf.lsferr.log',
    threads:
        config['tools']['bcftools']['threads']
    benchmark:
        BCFTOOLSOUT + '{sample}.vcf.benchmark'
    log:
        BCFTOOLSOUT + '{sample}.vcf.log'
    shell:
        ('{config[tools][bcftools][call]} ' +
        'call -A -v -c ' +
        '-o {output.vcf} ' +
        '{input.bcf} ' +
        '&& touch {output.suc}')

def getGATKVariantFiles():
    names = ''
    for w in SAMPLENAMES:
        names = names + '--variant ' + HAPLOTYPECALLEROUT + '/' + w + '.g.vcf '
    return names

if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = BASERECALIBRATIONOUT
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = OUTDIR + 'variants/GATK_HC/raw/'
rule gatkHaplotypeCaller:
    input:
        bam = HAPLOTYPECALLERIN + '{sample}.bam',
        bai = HAPLOTYPECALLERIN + '{sample}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = HAPLOTYPECALLEROUT + '{sample}.g.vcf',
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + '{sample}.g.vcf.lsferr.log',
        scratch = config['tools']['GATK']['haplotypeCaller']['scratch'],
        mem = config['tools']['GATK']['haplotypeCaller']['mem'],
        time = config['tools']['GATK']['haplotypeCaller']['time'],
    threads:
        config['tools']['GATK']['haplotypeCaller']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + '{sample}.g.vcf.log'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T HaplotypeCaller ' +
        '-R {input.reference} ' + 
        '-I {input.bam} ' +
        '--emitRefConfidence GVCF ' +
        '-L {input.regions} ' +
        '-o {output.vcf}')

rule gatkGenotypeGVCFs:
    input:
        vcf = expand(HAPLOTYPECALLEROUT + '{sample}.g.vcf', sample=SAMPLENAMES),
    output:
        vcf = HAPLOTYPECALLEROUT + 'combined.vcf',
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + 'combined.vcf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + 'combined.vcf.lsferr.log',
        scratch = config['tools']['GATK']['genotypeGVCFs']['scratch'],
        mem = config['tools']['GATK']['genotypeGVCFs']['mem'],
        time = config['tools']['GATK']['genotypeGVCFs']['time'],
        reference = config['resources'][ORGANISM]['reference'],
        input = getGATKVariantFiles()
    threads:
        config['tools']['GATK']['genotypeGVCFs']['threads']
    benchmark:
        HAPLOTYPECALLEROUT + 'combined.vcf.benchmark'
    log:
        HAPLOTYPECALLEROUT + 'combined.vcf.log'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T GenotypeGVCFs ' +
        '-R {params.reference} ' + 
        '{params.input} ' +
        '-o {output.vcf}')

# This is a GATK tool
# This is a GATK tool
if not 'MUTECT2IN' in globals():
    MUTECT2IN = BASERECALIBRATIONOUT
if not 'MUTECT2OUT' in globals():
    MUTECT2OUT = OUTDIR + 'variants/mutect2/raw/'
if not 'MUTECT2FILTEROUT' in globals():
    MUTECT2FILTEROUT = OUTDIR + 'variants/mutect2/filter/'
rule gatkMutect2:
    input:
        tumor = MUTECT2IN + '{tumor}.bam',
        tumorIdx = MUTECT2IN + '{tumor}.bai',
        normal = MUTECT2IN + '{normal}.bam',
        normalIdx = MUTECT2IN + '{normal}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        dbsnp = config['resources'][ORGANISM]['dbSNP'],
        cosmic = config['resources'][ORGANISM]['cosmic'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = MUTECT2OUT + '{tumor}_vs_{normal}.vcf',
        filterVcf = MUTECT2FILTEROUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = MUTECT2OUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = MUTECT2OUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['GATK']['mutect2']['scratch'],
        mem = config['tools']['GATK']['mutect2']['mem'],
        time = config['tools']['GATK']['mutect2']['time'],
        threads = config['tools']['GATK']['mutect2']['threads']
    threads:
        config['tools']['GATK']['mutect2']['threads']
    benchmark:
        MUTECT2OUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        MUTECT2OUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T MuTect2 ' +
        '-R {input.reference} ' + 
        '-I:tumor {input.tumor} ' +
        '-I:normal {input.normal} ' +
        '--dbsnp {input.dbsnp} ' +
        '--cosmic {input.cosmic} ' +
        '-L {input.regions} ' +
        '-nct {params.threads} ' +
        '-o {output.vcf}; ' + 
        'grep "^#" {output.vcf} > {output.filterVcf}; ' +
        'cat {output.vcf} | grep -v \"^#\" | awk \'{{if($7 == \"PASS\"){{print}}}}\' >> {output.filterVcf}')

def getDataBasisForMutect():
    out = []
    out.append(config['resources'][ORGANISM]['reference']) # this is a dummy such that something is retured
    if config['tools']['mutect1']['dbSNP']== "Y":
        out.append(config['resources'][ORGANISM]['dbSNP'])
    if config['tools']['mutect1']['cosmic']== "Y":
        out.append(config['resources'][ORGANISM]['cosmic'])
    return out

def prependDataBasisForRealignMutect():
    out = ""
    if config['tools']['mutect1']['dbSNP']== "Y":
        out += " --dbsnp " + config['resources'][ORGANISM]['dbSNP']
    if config['tools']['mutect1']['cosmic']== "Y":
        out += " --cosmic " + config['resources'][ORGANISM]['cosmic']
    return out

# Mutect1 is called separately from GATK
if not 'MUTECT1IN' in globals():
    MUTECT1IN = BASERECALIBRATIONOUT
if not 'MUTECT1OUT' in globals():
    MUTECT1OUT = OUTDIR + 'variants/mutect1/raw/'
rule mutect1:
    input:
        tumor = MUTECT1IN + '{tumor}.bam',
        tumorIdx = MUTECT1IN + '{tumor}.bai',
        normal = MUTECT1IN + '{normal}.bam',
        normalIdx = MUTECT1IN + '{normal}.bai',
        reference = {config['resources'][ORGANISM]['reference']},
        dbsnp = {config['resources'][ORGANISM]['dbSNP']},
        cosmic = {config['resources'][ORGANISM]['cosmic']},
        regions = {config['resources'][ORGANISM]['regions']}
    output:
        vcf = MUTECT1OUT + '{tumor}_vs_{normal}.vcf',
        out = MUTECT1OUT + '{tumor}_vs_{normal}.txt'
    params:
        lsfoutfile = MUTECT1OUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = MUTECT1OUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['mutect1']['scratch'],
        mem = config['tools']['mutect1']['mem'],
        time = config['tools']['mutect1']['time'],
        threads = config['tools']['mutect1']['threads'],
        intervalPadding = config['tools']['mutect1']['intervalPadding'],
        gapOpenPenalty = config['tools']['mutect1']['gapOpenPenalty'],
        tumorName = '{tumor}',
        normalName = '{normal}',
        dbs = prependDataBasisForRealignMutect()
    threads:
        config['tools']['mutect1']['threads']
    benchmark:
        MUTECT1OUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        MUTECT1OUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][mutect1][call]} ' +
        '--analysis_type MuTect ' +
        '--reference_sequence {input.reference} ' +
        '--input_file:tumor {input.tumor} ' +
        '--input_file:normal {input.normal} ' +
        '{params.dbs} ' +
        '--intervals {input.regions} ' +
        '--interval_padding {params.intervalPadding} ' +
        '--baqGapOpenPenalty {params.gapOpenPenalty} ' +
        '--out {output.out} ' +
        '--vcf {output.vcf}')

# call VarScan somatic
if not 'VARSCANSOMATICIN' in globals():
    VARSCANSOMATICIN = MPILEUPOUT
if not 'VARSCANSOMATICOUT' in globals():
    VARSCANSOMATICOUT = OUTDIR + 'variants/varscan_somatic/raw/'
rule varscanSomatic:
    input:
        tumor = VARSCANSOMATICIN + '{tumor}.mpileup',
        normal = VARSCANSOMATICIN + '{normal}.mpileup'
    output:
        vcfSnp = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.snp.vcf',
        vcfIndel = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.indel.vcf'
    params:
        lsfoutfile = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.lsfout.log',
        lsferrfile = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.lsferr.log',
        scratch = config['tools']['varscan']['somatic']['scratch'],
        mem = config['tools']['varscan']['somatic']['mem'],
        time = config['tools']['varscan']['somatic']['time'],
        outputTag = VARSCANSOMATICOUT + '{tumor}_vs_{normal}' 
    threads:
        config['tools']['varscan']['somatic']['threads']
    benchmark:
        VARSCANSOMATICOUT + '{tumor}_vs_{normal}.benchmark'
    log:
        VARSCANSOMATICOUT + '{tumor}_vs_{normal}.log'
    shell:
        ('{config[tools][varscan][call]} somatic {input.normal} {input.tumor} {params.outputTag} ' +
        '--output-vcf 1 ' +
        '{config[tools][varscan][somatic][params]}')
        
# call strelka, currently performed with the help of a shell script, needs to be changed soon!
# this script calls configureStrelkaWorkflow, then make on the temporary files. Then GATK is used to select only variants in the captured regions
# from snv and indel outputs. Finally, snv and indel output are combined
if not 'STRELKAIN' in globals():
    STRELKAIN = BASERECALIBRATIONOUT
if not 'STRELKAOUT' in globals():
    STRELKAOUT = OUTDIR + 'variants/strelka/'
rule strelka:
    input:
        tumor = STRELKAIN + '{tumor}.bam',
        tumorIdx = STRELKAIN + '{tumor}.bai',
        normal = STRELKAIN + '{normal}.bam',
        normalIdx = STRELKAIN + '{normal}.bai',
        dbSnpDB  = {config['resources'][ORGANISM]['dbSNP']},
        reference  = {config['resources'][ORGANISM]['reference']},
        regions = {config['resources'][ORGANISM]['regions']},
        strelkaConfig = {config['resources'][ORGANISM]['strelkaConfig']}
    output:
        #vcfSnp = STRELKAOUT + '{tumor}_vs_{normal}.snp.vcf',
        #vcfIndel = STRELKAOUT + '{tumor}_vs_{normal}.indel.vcf',
        vcf = STRELKAOUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = STRELKAOUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = STRELKAOUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['strelka']['scratch'],
        mem = config['tools']['strelka']['mem'],
        time = config['tools']['strelka']['time'],
        intervalPadding = config['tools']['strelka']['intervalPadding'],
        strelkaPerlScript = config['tools']['strelka']['strelkaPerlScript'],
        outDir = STRELKAOUT,
        outputTag = '{tumor}_vs_{normal}'
    threads:
        config['tools']['strelka']['threads']
    benchmark:
        STRELKAOUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        STRELKAOUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][strelka][call]} ' +
        '{params.outDir} ' +
        '{input.reference} ' + 
        '{input.dbSnpDB} ' +
        '{input.regions} ' +
        '{params.intervalPadding} ' +
        '{input.strelkaConfig} ' +
        '{input.tumor} ' +
        '{input.normal} ' +
        '{params.outputTag} ' +
        '{params.strelkaPerlScript} ' +
        '\"{config[tools][GATK][call]}\"')

# call somatic sniper
if not 'SOMATICSNIPERIN' in globals():
    SOMATICSNIPERIN = BASERECALIBRATIONOUT
if not 'SOMATICSNIPEROUT' in globals():
    SOMATICSNIPEROUT = OUTDIR + 'variants/somaticSniper/raw/'
rule somaticSniper:
    input:
        tumor = SOMATICSNIPERIN + '{tumor}.bam',
        tumorIdx = SOMATICSNIPERIN + '{tumor}.bai',
        normal = SOMATICSNIPERIN + '{normal}.bam',
        normalIdx = SOMATICSNIPERIN + '{normal}.bai',
        ref = config['resources'][ORGANISM]['reference']
    output:
        vcf = SOMATICSNIPEROUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = SOMATICSNIPEROUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = SOMATICSNIPEROUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['somaticSniper']['scratch'],
        mem = config['tools']['somaticSniper']['mem'],
        time = config['tools']['somaticSniper']['time'],
        outFormat = config['tools']['somaticSniper']['outFormat']
    threads:
        config['tools']['somaticSniper']['threads']
    benchmark:
        SOMATICSNIPEROUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        SOMATICSNIPEROUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][somaticSniper][call]} ' + 
        '-f {input.ref} ' +
        '{config[tools][somaticSniper][params]} ' + 
        '-F {params.outFormat} ' +
        '{input.tumor} ' +
        '{input.normal} ' +
        '{output.vcf}')

if not 'JOINTSNVMIX2_075_IN' in globals():
    JOINTSNVMIX2_075_IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2_075_OUT' in globals():
    JOINTSNVMIX2_075_OUT = OUTDIR + 'variants/jointSNVMix2_075/raw/'
rule JointSNVMix2_075_TRAIN:
    input:
        tumor = JOINTSNVMIX2_075_IN + '{tumor}.bam', 
        tumorIdx = JOINTSNVMIX2_075_IN + '{tumor}.bam.bai', 
        normal = JOINTSNVMIX2_075_IN + '{normal}.bam',
        normalIdx = JOINTSNVMIX2_075_IN + '{normal}.bam.bai',
        reference = {config['resources'][ORGANISM]['reference']},
        parameters = {config['resources']['general']['jsvm_0.7.5_jointParams']}, 
        priors = {config['resources']['general']['jsvm_0.7.5_jointPriors']}
    output:
        jsm = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_train.jsm'
    params:
        lsfoutfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_train.jsm.lsfout.log',
        lsferrfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_train.jsm.lsferr.log',
        scratch = config['tools']['jointSnvMix2_075']['train']['scratch'],
        mem = config['tools']['jointSnvMix2_075']['train']['mem'],
        time = config['tools']['jointSnvMix2_075']['train']['time'],
        skipSize = config['tools']['jointSnvMix2_075']['train']['skipSize'],
        convergence_threshold = config['tools']['jointSnvMix2_075']['train']['convergence_threshold'],
        method = config['tools']['jointSnvMix2_075']['method'],
    benchmark:
        JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_train.jsm.benchmark'
    threads:
        config['tools']['jointSnvMix2_075']['train']['threads']
    shell:
        ('{config[tools][jointSnvMix2_075][train][call]} ' +
        '{params.method} ' +
        '--skip_size {params.skipSize} ' +
        '{params.convergence_threshold} ' +
        '{input.reference} ' +
        '{input.normal} ' + 
        '{input.tumor} ' +
        '{input.priors} ' +
        '{input.parameters} ' +
        '{output.jsm}')

rule JointSNVMix2_075_CLASSIFY:
    input:
        tumor = JOINTSNVMIX2_075_IN + '{tumor}.bam', 
        tumorIdx = JOINTSNVMIX2_075_IN + '{tumor}.bam.bai', 
        normal = JOINTSNVMIX2_075_IN + '{normal}.bam',
        normalIdx = JOINTSNVMIX2_075_IN + '{normal}.bam.bai',
        reference = {config['resources'][ORGANISM]['reference']},
        parameters = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_train.jsm'
    output:
        jsm = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_classify.jsm'
    params:
        lsfoutfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_classify.jsm.lsfout.log',
        lsferrfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_classify.jsm.lsferr.log',
        scratch = config['tools']['jointSnvMix2_075']['classify']['scratch'],
        mem = config['tools']['jointSnvMix2_075']['classify']['mem'],
        time = config['tools']['jointSnvMix2_075']['classify']['time'],
        method = config['tools']['jointSnvMix2_075']['method']
    benchmark:
        JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_classify.jsm.benchmark'
    threads:
        config['tools']['jointSnvMix2_075']['classify']['threads']
    shell:
        ('{config[tools][jointSnvMix2_075][classify][call]} ' +
        '{params.method} ' +
        '{input.reference} ' +
        '{input.normal} ' +
        '{input.tumor} ' +
        '{input.parameters} ' +
        '{output.jsm}')

if not 'JOINTSNVMIX2IN' in globals():
    JOINTSNVMIX2IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2OUT' in globals():
    JOINTSNVMIX2OUT = OUTDIR + 'variants/jointSNVMix2/raw/'
rule JointSNVMix2_TRAIN:
    input:
        tumor = JOINTSNVMIX2IN + '{tumor}.bam', 
        tumorIdx = JOINTSNVMIX2IN + '{tumor}.bai', 
        normal = JOINTSNVMIX2IN + '{normal}.bam',
        normalIdx = JOINTSNVMIX2IN + '{normal}.bai',
        reference = {config['resources'][ORGANISM]['reference']},
        parameters = {config['resources']['general']['jsvmBetaBinParams']}, 
        priors = {config['resources']['general']['jsvmBetaBinPriors']}
    output:
        jsm = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_train.jsm'
    params:
        lsfoutfile = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_train.jsm.lsfout.log',
        lsferrfile = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_train.jsm.lsferr.log',
        scratch = config['tools']['jointSnvMix2']['train']['scratch'],
        mem = config['tools']['jointSnvMix2']['train']['mem'],
        time = config['tools']['jointSnvMix2']['train']['time'],
        #minBaseQual = '--min_base_qual ' + config['tools']['jointSnvMix2']['minBaseQual'],
        minBaseQual = '',
        skipSize = config['tools']['jointSnvMix2']['train']['skipSize'],
        params = config['tools']['jointSnvMix2']['train']['params']
    benchmark:
        JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_train.jsm.benchmark'
    threads:
        config['tools']['jointSnvMix2']['train']['threads']
    shell:
        ('{config[tools][jointSnvMix2][train][call]} ' +
        '{params.minBaseQual} ' +
        '--skip_size {params.skipSize} ' +
        '{params.params} ' +
        '--model beta_binomial ' +
        '--priors {input.priors} ' +
        '--initial_parameters_file {input.parameters} ' +
        '{input.reference} ' +
        '{input.normal} ' + 
        '{input.tumor} ' +
        '{output.jsm}')

rule JointSNVMix2_CLASSIFY:
    input:
        tumor = JOINTSNVMIX2IN + '{tumor}.bam', 
        normal = JOINTSNVMIX2IN + '{normal}.bam',
        reference = {config['resources'][ORGANISM]['reference']},
        parameters = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_train.jsm'
    output:
        jsm = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_classify.jsm',
        vcf = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_classify.jsm.lsfout.log',
        lsferrfile = JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_classify.jsm.lsferr.log',
        scratch = config['tools']['jointSnvMix2']['classify']['scratch'],
        mem = config['tools']['jointSnvMix2']['classify']['mem'],
        time = config['tools']['jointSnvMix2']['classify']['time'],
        minBaseQual = ''
        #minBaseQual = '--min_base_qual ' + config['tools']['jointSnvMix2']['minBaseQual']
    benchmark:
        JOINTSNVMIX2OUT + '{tumor}_vs_{normal}_classify.jsm.benchmark'
    threads:
        config['tools']['jointSnvMix2']['classify']['threads']
    shell:
        ('{config[tools][jointSnvMix2][classify][call]} ' +
        '{params.minBaseQual} ' +
        '--model beta_binomial ' +
        '--out_file {output.jsm} ' +
        '--parameters_file {input.parameters} ' +
        '--post_process ' +
        '{input.reference} ' +
        '{input.normal} ' +
        '{input.tumor}; '
        '{config[tools][jointSnvMix2][jsm2Vcf][call]} {output.jsm} {output.vcf} {config[tools][jointSnvMix2][jsm2Vcf][minProb]}')

if not 'VARDICTIN' in globals():
    VARDICTIN = BASERECALIBRATIONOUT
if not 'VARDICTOUT' in globals():
    VARDICTOUT = OUTDIR + 'variants/varDict/raw/'
if not 'VARDICTFILTEROUT' in globals():
    VARDICTFILTEROUT = OUTDIR + 'variants/varDict/filtered/'
rule varDict:
    input:
        tumor = VARDICTIN + '{tumor}.bam',
        tumorIdx = VARDICTIN + '{tumor}.bai',
        normal = VARDICTIN + '{normal}.bam',
        normalIdx = VARDICTIN + '{normal}.bai',
        regions = config['resources'][ORGANISM]['regions'],
        ref = config['resources'][ORGANISM]['reference']
    output:
        vcf = VARDICTOUT + '{tumor}_vs_{normal}.vcf',
        filterVcf = VARDICTFILTEROUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = VARDICTOUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = VARDICTOUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['varDict']['scratch'],
        mem = config['tools']['varDict']['mem'],
        time = config['tools']['varDict']['time'],
        minAllelFreq = config['tools']['varDict']['minAllelFreq'],
        tumor = '{tumor}',
        normal = '{normal}'
    threads:
        config['tools']['varDict']['threads']
    benchmark:
        VARDICTOUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        VARDICTOUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][varDict][expRscript]} ' + 
        '{config[tools][varDict][call]} ' + 
        '-G {input.ref} ' +
        '-f {params.minAllelFreq} ' +
        '-b "{input.tumor}|{input.normal}" ' +
        '-Q 1 -c 1 -S 2 -E 3 -g 4 {input.regions} | ' +
        '{config[tools][varDict][varDictTestSomatic]} | ' +
        '{config[tools][varDict][varDict2VcfPaired]} ' +
        '-N "TUMOR|NORMAL" ' +
        '-f {params.minAllelFreq} '
        '> {output.vcf}; ' +
        'grep "^#" {output.vcf} > {output.filterVcf}; ' +
        'cat {output.vcf} | grep -v \"^#\" | awk \'{{if($7 == \"PASS\"){{print}}}}\' >> {output.filterVcf}')


if not 'DEEPSNVIN' in globals():
    DEEPSNVIN = BASERECALIBRATIONOUT
if not 'DEEPSNVOUT' in globals():
    DEEPSNVOUT = OUTDIR + 'variants/deepSNV/raw/'
rule deepSNV:
    input:
        tumor = DEEPSNVIN + '{tumor}.bam',
        tumorIdx = DEEPSNVIN + '{tumor}.bai',
        normal = DEEPSNVIN + '{normal}.bam',
        normalIdx = DEEPSNVIN + '{normal}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = DEEPSNVOUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = DEEPSNVOUT + '{tumor}_vs_{normal}.vcf.lsfout.log',
        lsferrfile = DEEPSNVOUT + '{tumor}_vs_{normal}.vcf.lsferr.log',
        scratch = config['tools']['deepSNV']['scratch'],
        mem = config['tools']['deepSNV']['mem'],
        time = config['tools']['deepSNV']['time'],
        outDir = DEEPSNVOUT,
        params = config['tools']['deepSNV']['params']
    threads:
        config['tools']['deepSNV']['threads']
    benchmark:
        DEEPSNVOUT + '{tumor}_vs_{normal}.vcf.benchmark'
    log:
        DEEPSNVOUT + '{tumor}_vs_{normal}.vcf.log'
    shell:
        ('{config[tools][deepSNV][call]} ' + 
        '{input.regions} ' +
        '{input.tumor} ' +
        '{input.normal} ' +
        '{input.reference} ' +
        '{params.outDir} ' +
        '{params.params} ' +
        '--also-vcf true')

def getVcfsForRankComination(wildcards):
    out = []
    if config['tools']['rankCombineVariants']['mutect'] == "Y":
        out.append(MUTECT1OUT + wildcards.sample +'.vcf')
    if config['tools']['rankCombineVariants']['deepSNV'] == "Y":
        out.append(DEEPSNVOUT + wildcards.sample +'.vcf')
    if config['tools']['rankCombineVariants']['varscan'] == "Y":
        out.append(VARSCANSOMATICOUT + wildcards.sample +'.vcf')
    if config['tools']['rankCombineVariants']['jsvm'] == "Y":
        out.append(JOINTSNVMIX2OUT + wildcards.sample +'.vcf')
    if config['tools']['rankCombineVariants']['strelka'] == "Y":
        out.append(STRELKAOUT + wildcards.sample +'.vcf')
    return out

if not 'RANKCOMBINEOUT' in globals():
    RANKCOMBINEOUT = OUTDIR + 'variants/rank_combine/'
rule rankCombineVariants:
    input:
        vcfs = getVcfsForRankComination
    output:
        txt = RANKCOMBINEOUT + '{sample}.txt'
    params:
        lsfoutfile = RANKCOMBINEOUT + '{sample}.vcf.lsfout.log',
        lsferrfile = RANKCOMBINEOUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['rankCombineVariants']['scratch'],
        mem = config['tools']['rankCombineVariants']['mem'],
        time = config['tools']['rankCombineVariants']['time']
    threads:
        config['tools']['rankCombineVariants']['threads']
    benchmark:
        RANKCOMBINEOUT + '{sample}.vcf.benchmark'
    log:
        RANKCOMBINEOUT + '{sample}.vcf.log'
    shell:
        ('{config[tools][rankCombineVariants][call]} ' + 
        '{output.txt} ' + 
        '{input.vcfs}')

