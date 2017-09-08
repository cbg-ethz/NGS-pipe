from snakemake.utils import R

def getContigNames():
    if isinstance(config['resources'][ORGANISM]['contigNames'], Error):
        return ["ERROR"]
    f = open(config['resources'][ORGANISM]['contigNames'],'r')
    result = []
    for line in f:
        result.append(line.strip())
    return result

# This rule uses the custom samtools version of bicSeq2 to extract unique mappings for bicSeq2
if not 'BICSEQ2IN' in globals():
    BICSEQ2IN = REMOVEPCRDUBLICATESOUT
if not 'BICSEQ2OUT' in globals():
    BICSEQ2OUT = OUTDIR + 'bicseq2/'
rule bicSeq_samtoolsUnique:
    input:
        bam = BICSEQ2IN + '{sample}.bam',
        contigNamnes =config['resources'][ORGANISM]['contigNames']
    output:
        seq = expand(BICSEQ2OUT + '{{sample}}' + '/{contigNames}.seq', contigNames=getContigNames())
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}/unique.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}/unique.lsferr.log',
        scratch = config['tools']['bicseq2']['unique']['scratch'],
        mem = config['tools']['bicseq2']['unique']['mem'],
        time = config['tools']['bicseq2']['unique']['time'],
        mapper = config['tools']['bicseq2']['unique']['mapper'],
        directory = BICSEQ2OUT + '{sample}' + '/'
    threads:
        config['tools']['bicseq2']['unique']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}/unique.benchmark'
    shell:
        ('mkdir -p {params.directory} && ' +
        '{config[tools][bicseq2][unique][call]} view -U {params.mapper},{params.directory},N,N {input.bam}')

rule extractContigsFromFasta:
    input:
        fasta = config['resources'][ORGANISM]['reference']
    output:
        fasta = expand(config['resources'][ORGANISM]['reference'] + '_contigs/{contigs}.fasta', contigs = getContigNames()),
        suc = config['resources'][ORGANISM]['reference'] + '_contigs/complete.txt'
    params:
        lsfoutfile = config['resources'][ORGANISM]['reference'] + '_contigs/complete.lsfout.log',
        lsferrfile = config['resources'][ORGANISM]['reference'] + '_contigs/complete.lsferr.log',
        scratch = config['tools']['extractContigs']['scratch'],
        mem = config['tools']['extractContigs']['mem'],
        time = config['tools']['extractContigs']['time'],
    threads:
        config['tools']['extractContigs']['threads']
    benchmark:
        config['resources'][ORGANISM]['reference'] + '_contigs/complete.benchmark'
    shell:
        '{config[tools][extractContigs][call]} {input.fasta} && touch {output.suc}'


rule craeteConfigBicSeqNorm:
    input:
        suc = config['resources'][ORGANISM]['reference'] + '_contigs/complete.txt'
    output:
        config = BICSEQ2OUT + '{sample}/configNorm.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}/configNorm.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}/configNorm.txt.lsferr.log',
        scratch = config['tools']['bicSeqConfigNorm']['scratch'],
        mem = config['tools']['bicSeqConfigNorm']['mem'],
        time = config['tools']['bicSeqConfigNorm']['time'],
        sample = '{sample}',
        reference = config['resources'][ORGANISM]['reference'],
        mappabilityFile = config['resources'][ORGANISM]['pathBicSeq2Mappability']
    threads:
        config['tools']['bicSeqConfigNorm']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}/configNorm.txt.benchmark'
    run:
        outfile = open(output.config, 'w')
        outfile.write('chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n')
        #TODO: This needs to be adapted to work for different organisms and versions
        BICSEQMAPPINGCHROM = sorted([file.replace(config['resources'][ORGANISM]['pathBicSeq2Mappability'], '').strip().split('.')[-2] for file in glob.glob(config['resources'][ORGANISM]['pathBicSeq2Mappability'] + '/*.txt')])
        for chr in BICSEQMAPPINGCHROM:
            outfile.write(chr + '\t')
            outfile.write(params.reference + '_contigs/' + chr + '.fasta' + '\t')
            outfile.write(params.mappabilityFile + '/hg19.CRC.75mer.' + chr + '.txt' + '\t')
            outfile.write(BICSEQ2OUT + params.sample + '/' + chr + '.seq\t')
            outfile.write(BICSEQ2OUT + params.sample + '/' + chr + '.norm.bin\n')

# This rule applies bicSeq2-norm
rule bicSeq_norm:
    input:
        insertSizeFile = BICSEQ2IN + '{sample}.bam_stats/insert_size_metrics.txt',
        config = BICSEQ2OUT + '{sample}/configNorm.txt',
        seq = expand(BICSEQ2OUT + '{{sample}}' + '/{contigNames}.seq', contigNames=getContigNames())
    output:
        out = BICSEQ2OUT + '{sample}/paramsEstimate.txt',
        tmp = temp(BICSEQ2OUT + '{sample}/tmp/')
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}/paramsEstimate.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}/paramsEstimate.txt.lsferr.log',
        scratch = config['tools']['bicseq2']['norm']['scratch'],
        mem = config['tools']['bicseq2']['norm']['mem'],
        time = config['tools']['bicseq2']['norm']['time'],
        readLength = config['tools']['bicseq2']['norm']['readLength'],
        directory = BICSEQ2OUT + '{sample}'
    threads:
        config['tools']['bicseq2']['norm']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}/paramsEstimate.txt.benchmark'
    shell:
        ("inSize=$(head -n 8 {input.insertSizeFile} | tail -n 1 | cut -f 1) && " +
        "{config[tools][bicseq2][norm][call]} " +
        "-l={params.readLength} " +
        "-s=${{inSize}} " +
        "--gc_bin {input.config} " +
        "--tmp={output.tmp} " +
        "{output.out}")

rule craeteConfigBicSeqSeg:
    input:
        suc = config['resources'][ORGANISM]['reference'] + '_contigs/complete.txt'
    output:
        config = BICSEQ2OUT + '{tumor}_vs_{normal}/configSeg.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{tumor}_vs_{normal}/configSeg.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{tumor}_vs_{normal}/configSeg.txt.lsferr.log',
        scratch = config['tools']['bicSeqConfigSeq']['scratch'],
        mem = config['tools']['bicSeqConfigSeq']['mem'],
        time = config['tools']['bicSeqConfigSeq']['time'],
        normal = '{normal}',
        tumor = '{tumor}'
    threads:
        config['tools']['bicSeqConfigSeq']['threads']
    benchmark:
        BICSEQ2OUT + '{tumor}_vs_{normal}/configSeg.txt.benchmark'
    run:
        outfile = open(output.config, 'w')
        outfile.write('chromName\tbinFileNorm.Case\tbinFileNorm.Control\n')
        #TODO: This needs to be adapted to work for different organisms and versions
        BICSEQMAPPINGCHROM = sorted([file.replace(config['resources'][ORGANISM]['pathBicSeq2Mappability'], '').strip().split('.')[-2] for file in glob.glob(config['resources'][ORGANISM]['pathBicSeq2Mappability'] + '/*.txt')])
        for chr in BICSEQMAPPINGCHROM:
            outfile.write(chr + '\t')
            outfile.write(BICSEQ2OUT + params.tumor + '/' + chr + '.norm.bin\t')
            outfile.write(BICSEQ2OUT + params.normal + '/' + chr + '.norm.bin\n')

# This rule applies bicSeq2-seg
rule bicSeq_seg:
    input:
        config = BICSEQ2OUT + '{tumor}_vs_{normal}/configSeg.txt',
        paramsNormal = BICSEQ2OUT + '{normal}/paramsEstimate.txt',
        paramsTumor = BICSEQ2OUT + '{tumor}/paramsEstimate.txt'
    output:
        out = BICSEQ2OUT + '{tumor}_vs_{normal}.cnvsRaw.txt',
        fig = BICSEQ2OUT + '{tumor}_vs_{normal}.png',
        tmp = temp(BICSEQ2OUT + '{tumor}_vs_{normal}/tmp/')
    params:
        lsfoutfile = BICSEQ2OUT + '{tumor}_vs_{normal}.cnvsRaw.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{tumor}_vs_{normal}.cnvsRaw.txt.lsferr.log',
        scratch = config['tools']['bicseq2']['seg']['scratch'],
        mem = config['tools']['bicseq2']['seg']['mem'],
        time = config['tools']['bicseq2']['seg']['time'],
        title = '{tumor}_vs_{normal}_CNV'
    threads:
        config['tools']['bicseq2']['seg']['threads']
    benchmark:
        BICSEQ2OUT + '{tumor}_vs_{normal}.cnvsRaw.txt.benchmark'
    shell:
        '{config[tools][bicseq2][seg][call]} ' +
        '--fig={output.fig} ' +
        '--title={params.title} ' +
        '--nrm ' + 
        '--control ' +
        '--tmp={output.tmp} ' +
        '{input.config} {output.out}'
        #'{config[tools][bicseq2][seg][call]} --fig={output.fig} --title={params.title} --control {input.config} {output.out}'


# This rule applies bicSeq2 genotype.pl to assess event significance
rule bicSeq_genotype:
    input:
        inCNV = BICSEQ2OUT + '{sample}.cnvsRaw.txt',
        config = BICSEQ2OUT + '{sample}/configSeg.txt'
    output:
        out = BICSEQ2OUT + '{sample}.cnvsGenotype.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}.cnvsGenotype.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}.cnvsGenotype.txt.lsferr.log',
        scratch = config['tools']['bicseq2']['genotype']['scratch'],
        mem = config['tools']['bicseq2']['genotype']['mem'],
        time = config['tools']['bicseq2']['genotype']['time']
    threads:
        config['tools']['bicseq2']['genotype']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}.cnvsGenotype.txt.benchmark'
    shell:
        'cut -f 1-3 {input.inCNV} > {input.inCNV}.forGenotype && ' +
        '{config[tools][bicseq2][genotype][call]} {input.config} {input.inCNV}.forGenotype {output.out}'


# This rule applies a simple filter script to filter cnv events given a certain pvalue threshold and to determine the copy number
rule bicSeq_filter:
    input:
        inCNV = BICSEQ2OUT + '{sample}.cnvsGenotype.txt'
    output:
        out = BICSEQ2OUT + '{sample}.filtered.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}.filtered.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}.filtered.txt.lsferr.log',
        scratch = config['tools']['bicseq2']['filter']['scratch'],
        mem = config['tools']['bicseq2']['filter']['mem'],
        time = config['tools']['bicseq2']['filter']['time'],
        pvalue = config['tools']['bicseq2']['filter']['pvalueThreshold']
    threads:
        config['tools']['bicseq2']['filter']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}.filtered.txt.benchmark'
    shell:
        '{config[tools][bicseq2][filter][call]} {input.inCNV} {output.out} {params.pvalue}'


# call VarScan copynumner
if not 'VARSCANCNVIN' in globals():
    VARSCANCNVIN = MPILEUPOUT
if not 'VARSCANCNVOUT' in globals():
    VARSCANCNVOUT = OUTDIR + 'copy_number/varscan_cnv/'
rule varscan_copy_number:
    input:
        tumor = VARSCANCNVIN + '{tumor}.mpileup',
        normal = VARSCANCNVIN + '{normal}.mpileup'
    output:
        out = VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber'
    params:
        lsfoutfile = VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber.lsfout.log',
        lsferrfile = VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber.lsferr.log',
        scratch = config['tools']['varscan']['copyNumber']['scratch'],
        mem = config['tools']['varscan']['copyNumber']['mem'],
        time = config['tools']['varscan']['copyNumber']['time'],
        params = config['tools']['varscan']['copyNumber']['params'],
        outputTag = VARSCANCNVOUT + '{tumor}_vs_{normal}' 
    threads:
        config['tools']['varscan']['copyNumber']['threads']
    benchmark:
        VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber.benchmark'
    log:
        VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber.log'
    shell:
        ('{config[tools][varscan][call]} copynumber ' +
        '{input.normal} ' +
        '{input.tumor} ' +
        '{params.outputTag} ' +
        '{params.params}')

# call VarScan copyCaller
rule varscan_copy_caller:
    input:
        rawCN = VARSCANCNVOUT + '{tumor}_vs_{normal}.copynumber'
    output:
        out = VARSCANCNVOUT + '{tumor}_vs_{normal}.cn'
    params:
        lsfoutfile = VARSCANCNVOUT + '{tumor}_vs_{normal}.cn.lsfout.log',
        lsferrfile = VARSCANCNVOUT + '{tumor}_vs_{normal}.cn.lsferr.log',
        scratch = config['tools']['varscan']['copyCaller']['scratch'],
        mem = config['tools']['varscan']['copyCaller']['mem'],
        time = config['tools']['varscan']['copyCaller']['time'],
        params = config['tools']['varscan']['copyCaller']['params']
    threads:
        config['tools']['varscan']['copyCaller']['threads']
    benchmark:
        VARSCANCNVOUT + '{tumor}_vs_{normal}.cy.benchmark'
    log:
        VARSCANCNVOUT + '{tumor}_vs_{normal}.cn.log'
    shell:
        ('{config[tools][varscan][call]} copyCaller ' +
        '{input.rawCN} ' + 
        '--output-file {output.out} ' +
        '{params.params}')

rule bicSeq2annovar:
    input:
        BICSEQ2OUT + '{sample}.filtered.txt'
    output:
        BICSEQ2OUT + '{sample}.filtered.forAnnovar.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}.filtered.forAnnovar.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}.filtered.forAnnovar.txt.lsferr.log',
        scratch = config['tools']['bicSeq2annovar']['scratch'],
        mem = config['tools']['bicSeq2annovar']['mem'],
        time = config['tools']['bicSeq2annovar']['time']
    threads:
        config['tools']['bicSeq2annovar']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}.filtered.forAnnovar.txt.benchmark'
    shell:
        '{config[tools][bicSeq2annovar][call]} {input} {output}'

#TODO: This needs to be adapted to work for different organisms and versions
rule wgsAnnovar:
    input:
        txt = BICSEQ2OUT + '{sample}.filtered.forAnnovar.txt',
        db = config['resources'][ORGANISM]['annovarDB']
    output:
        out = BICSEQ2OUT + '{sample}.filtered.annotated.hg19_multianno.txt'
    params:
        lsfoutfile = BICSEQ2OUT + '{sample}.filtered.annotated.hg19_multianno.txt.lsfout.log',
        lsferrfile = BICSEQ2OUT + '{sample}.filtered.annotated.hg19_multianno.txt.lsferr.log',
        scratch = config['tools']['annovar']['scratch'],
        mem = config['tools']['annovar']['mem'],
        time = config['tools']['annovar']['time'],
        buildver = config['tools']['annovar']['buildver'],
        params = config['tools']['annovar']['params'],
        out = BICSEQ2OUT + '{sample}.filtered.annotated'
    threads:
        config['tools']['annovar']['threads']
    benchmark:
        BICSEQ2OUT + '{sample}.filtered.annotated.hg19_multianno.txt.benchmark'
    shell:
        ('{config[tools][annovar][call]} ' +
        '{input.txt} ' +
        '{input.db} ' +
        '-buildver {params.buildver} ' +
        '-out {params.out} ' +
        '{params.params}')

if not 'FACETSIN' in globals():
    FACETSIN = BASERECALIBRATIONOUT
if not 'FACETSOUT' in globals():
    FACETSOUT = OUTDIR + 'copy_number/facets/'
rule createBedForFacets:
    input:
        vcf = config['resources'][ORGANISM]['dbSNP'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        vcf = FACETSOUT + 'snps.vcf'
    params:
        lsfoutfile = FACETSOUT + 'snps.vcf.lsfout.log',
        lsferrfile = FACETSOUT + 'snps.vcf.lsferr.log',
        scratch = config['tools']['facets']['region']['scratch'],
        mem = config['tools']['facets']['region']['mem'],
        time = config['tools']['facets']['region']['time'],
        params = config['tools']['facets']['region']['params']
    threads:
        config['tools']['facets']['region']['threads']
    benchmark:
        FACETSOUT + 'snps.vcf.benchmark'
    shell:
        ('grep "^#" {input.vcf} > {output.vcf}; ' +
        '{config[tools][facets][region][call]} ' +
        '{params.params} ' +
        '-a {input.vcf} ' +
        '-b {input.regions} ' +
        '>> {output.vcf}')

rule getSNPInfoForFacets:
    input:
        vcf = FACETSOUT + 'snps.vcf',
        normal = FACETSIN + '{normal}.bam',
        tumor = FACETSIN + '{tumor}.bam'
    output:
        csv = FACETSOUT + '{tumor}_vs_{normal}.csv.gz'
    params:
        lsfoutfile = FACETSOUT + '{tumor}_vs_{normal}.csv.gz.lsfout.log',
        lsferrfile = FACETSOUT + '{tumor}_vs_{normal}.csv.gz.lsferr.log',
        scratch = config['tools']['facets']['snpPileup']['scratch'],
        mem = config['tools']['facets']['snpPileup']['mem'],
        time = config['tools']['facets']['snpPileup']['time'],
        params = config['tools']['facets']['snpPileup']['params']
    threads:
        config['tools']['facets']['snpPileup']['threads']
    benchmark:
        FACETSOUT + '{tumor}_vs_{normal}.csv.gz.benchmark'
    shell:
        ('{config[tools][facets][snpPileup][call]} ' +
        '{params.params} ' +
        '{input.vcf} ' +
        '{output.csv} ' +
        '{input.normal} ' +
        '{input.tumor}')

rule facets:
    input:
        csv = FACETSOUT + '{tumor}_vs_{normal}.csv.gz'
    output:
        pdf = FACETSOUT + '{tumor}_vs_{normal}.pdf',
        txt = FACETSOUT + '{tumor}_vs_{normal}.cn'
    params:
        lsfoutfile = FACETSOUT + '{tumor}_vs_{normal}.cn.lsfout.log',
        lsferrfile = FACETSOUT + '{tumor}_vs_{normal}.cn.lsferr.log',
        scratch = config['tools']['facets']['facets']['scratch'],
        mem = config['tools']['facets']['facets']['mem'],
        time = config['tools']['facets']['facets']['time'],
        params = config['tools']['facets']['facets']['params']
    threads:
        config['tools']['facets']['facets']['threads']
    benchmark:
        FACETSOUT + '{tumor}_vs_{normal}.cn.benchmark'
    shell:
        ('{config[tools][facets][facets][call]} ' + 
        '{params.params} ' +
        '{input.csv} {output.txt} {output.pdf}')

# This rule annotates the CNV call results (for excavator reformatting is needed first)
# TODO: may be of general use, if database is part of paramters? Chose a better suited place in that case!
rule annotateCNVsWithBedtools:
    input:
        inRes = '{sample}.txt',
        inDB = config['resources']['H_sapiens_hg19']['geneAnnotationDB'] 
    output:
        out = '{sample}.annotated.txt'
    params:
        lsfoutfile = '{sample}.annotated.lsfout.log',
        lsferrfile = '{sample}.annotated.lsferr.log',
        scratch = config['tools']['bedtools']['intersect']['scratch'],
        mem = config['tools']['bedtools']['intersect']['mem'],
        time = config['tools']['bedtools']['intersect']['time']
    threads:
        config['tools']['bedtools']['intersect']['threads']
    benchmark:
        '{sample}.annotated.benchmark'
    shell:
        '{config[tools][bedtools][intersect][call]} intersect -a {input.inRes} -b {input.inDB} -wa -wb > {output.out}.temp1.txt ; ' +
        '{config[tools][bedtools][intersect][call]} intersect -a {input.inRes} -b {input.inDB} -v -wa > {output.out}.temp2.txt ; ' +
        'cat {output.out}.temp1.txt {output.out}.temp2.txt | sort -k1,1 -k2,2n > {output.out}'

