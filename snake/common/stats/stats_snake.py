# This rule uses qualimap to generate a report about several statistics of a BAM file.
localrules: createBedForQualimap
rule createBedForQualimap:
    input:
        regions = config['resources'][ORGANISM]['regions']
    output:
        regions = config['resources'][ORGANISM]['regions'] + '_qual.bed'
    shell:
        'awk \'{{sub(/\r/,\"\"); if(NF > 1){{printf $0; for(i = NF; i < 6; ++i){{printf \"\\t*\"}}; printf \"\\n\"}}}}\' {input.regions} > {output.regions}'

rule qualimap_PDF:
    input:
        bam = '{sample}.bam',
        regions = config['resources'][ORGANISM]['regions'] + '_qual.bed' 
    output:
        dir = '{sample}.bam_stats',
        file = '{sample}.bam_stats/report.pdf'
    params:
        lsfoutfile = '{sample}.bam_stats/report.pdf.lsfout.log',
        lsferrfile = '{sample}.bam_stats/report.pdf.lsferr.log',
        scratch = config['tools']['qualimap']['scratch'],
        mem = config['tools']['qualimap']['mem'],
        time = config['tools']['qualimap']['time'],
        params = config['tools']['qualimap']['params']
    benchmark:
        '{sample}.bam_stats/report.pdf.benchmark'
    threads:
        config['tools']['qualimap']['threads']
    shell:
        'if [[ ! -n $({config[tools][samtools][call]} view {input.bam} | head -n 1) ]]; then touch {output.file}; else {config[tools][qualimap][call]} bamqc -bam {input.bam} -outdir {output.dir} -outformat PDF -os -feature-file {input.regions} --java-mem-size={config[tools][qualimap][mem]}M {params.params}; fi'

rule qualimap_HTML:
    input:
        bam = '{sample}.bam',
        regions = config['resources'][ORGANISM]['regions'] + '_qual.bed'
    output:
        dir = '{sample}.bam_stats',
        file = '{sample}.bam_stats/qualimapReport.html'
    params:
        lsfoutfile = '{sample}.bam_stats/qualimapReport.html.lsfout.log',
        lsferrfile = '{sample}.bam_stats/qualimapReport.html.lsferr.log',
        scratch = config['tools']['qualimap']['scratch'],
        mem = config['tools']['qualimap']['mem'],
        time = config['tools']['qualimap']['time'],
        params = config['tools']['qualimap']['params']
    benchmark:
        '{sample}.bam_stats/qualimapReport.html.benchmark'
    threads:
        config['tools']['qualimap']['threads']
    shell:
        'if [[ ! -n $({config[tools][samtools][call]} view {input.bam} | head -n 1) ]]; then touch {output.file}; else {config[tools][qualimap][call]} bamqc -bam {input.bam} -outdir {output.dir} -outformat HTML -os -feature-file {input.regions} --java-mem-size={config[tools][qualimap][mem]}M {params.params}; fi'

rule multiqcBam:
    input:
        qualimap = expand('{{analysisStep}}/{sample}.bam_stats/qualimapReport.html', sample = SAMPLENAMES),
        samtools = expand('{{analysisStep}}/{sample}.bam.flagstat', sample = SAMPLENAMES),
    output:
        '{analysisStep}/multiqc_report.html'
    params:
        lsfoutfile = '{analysisStep}/multiqc_report.html.lsfout.log',
        lsferrfile = '{analysisStep}/multiqc_report.html.lsferr.log',
        scratch = config['tools']['multiqc']['scratch'],
        mem = config['tools']['multiqc']['mem'],
        time = config['tools']['multiqc']['time'],
        params = config['tools']['multiqc']['params']
    benchmark:
        '{analysisStep}/multiqc_report.html.benchmark'
    threads:
        config['tools']['multiqc']['threads']
    shell:
        '{config[tools][multiqc][call]} --outdir {wildcards.analysisStep} --filename multiqc_report.html {params.params} {wildcards.analysisStep}'

rule multiqcFastq:
    input:
        fastqc = expand('{{analysisStep}}/{fastq}_fastqc.html', fastq = getPairedFastqFiles(SAMPLENAMES)),
    output:
        '{analysisStep}/multiqc_report.html'
    params:
        lsfoutfile = '{analysisStep}/multiqc_report.html.lsfout.log',
        lsferrfile = '{analysisStep}/multiqc_report.html.lsferr.log',
        scratch = config['tools']['multiqc']['scratch'],
        mem = config['tools']['multiqc']['mem'],
        time = config['tools']['multiqc']['time'],
        params = config['tools']['multiqc']['params']
    benchmark:
        '{analysisStep}/multiqc_report.html.benchmark'
    threads:
        config['tools']['multiqc']['threads']
    shell:
        '{config[tools][multiqc][call]} --ignore *R3_fastqc.html --outdir {wildcards.analysisStep} --filename multiqc_report.html {params.params} {wildcards.analysisStep}'

def exportRScript(path):
    if path != "":
        return path + ";"
    else:
        return ""

# This rule gives statistics about the effect of GATKs base recalibration had on the BAM file.
rule analyzeCovariates:
    input:
        tabBefore = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
        tabAfter = BASERECALIBRATIONOUT + '{sample}.secondPass_reca.table',
        reference = config['resources'][ORGANISM]['reference']
    output:
        pdf=BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf'
    params:
        lsfoutfile = BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf.lsfout.log',
        lsferrfile = BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf.lsferr.log',
        scratch = config['tools']['GATK']['analyzeCovariates']['scratch'],
        mem = config['tools']['GATK']['analyzeCovariates']['mem'],
        time = config['tools']['GATK']['analyzeCovariates']['time'],
        params = config['tools']['GATK']['analyzeCovariates']['params'],
        rScriptPath = exportRScript(config['tools']['GATK']['analyzeCovariates']['expRscript'])
    benchmark:
        BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf.benchmark'
    threads:
        config['tools']['GATK']['analyzeCovariates']['threads']
    shell:
        ('{params.rScriptPath} ' +
        '{config[tools][GATK][call]} ' +
        '-T AnalyzeCovariates ' +
        '-R {input.reference} ' +
        '-before {input.tabBefore} ' +
        '-after {input.tabAfter} ' +
        '-plots {output.pdf} ' +
        '{params.params}')

rule countReadsInFastq:
    input:
        fastq = '{sample}.fastq.gz'
    output:
        count = '{sample}.count'
    params:
        mem = '1000',
        scratch = '1000',
        time = '20',
        lsfoutfile = '{sample}.count.lsfout.log',
        lsferrfile = '{sample}.count.lsferr.log'
    benchmark:
        '{sample}.count.benchmark'
    threads:
        1
    shell:
        'zcat {input.fastq} | wc -l | awk \'{{print $1/4}}\' > {output.count}'

rule samtools_flagstat:
    input:
        bam = '{sample}.bam',
    output:
        flagstat = '{sample}.bam.flagstat',
    params:
        mem = config['tools']['samtools']['flagstat']['mem'],
        scratch = config['tools']['samtools']['flagstat']['scratch'],
        time = config['tools']['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.bam.flagstat.lsfout.log',
        lsferrfile = '{sample}.bam.flagstat.lsferr.log',
        params = config['tools']['samtools']['flagstat']['params']
    threads:
        config['tools']['samtools']['flagstat']['threads']
    shell:
        '{config[tools][samtools][call]} flagstat {params.params} {input.bam} > {output.flagstat}'

# This rule uses qualimap to generate a report about several statistics of a BAM file.
rule insertSizeMetric:
    input:
        bam = '{sample}.bam'
    output:
        dir = '{sample}.bam_stats',
        txt = '{sample}.bam_stats/insert_size_metrics.txt',
        pdf = '{sample}.bam_stats/insert_size_histogram.pdf'
    params:
        lsfoutfile = '{sample}.bam_stats/insert_size_metrics.txt.lsfout.log',
        lsferrfile = '{sample}.bam_stats/insert_size_metrics.txt.lsferr.log',
        scratch = config['tools']['picard']['collectInsertSizeMetrics']['scratch'],
        mem = config['tools']['picard']['collectInsertSizeMetrics']['mem'],
        time = config['tools']['picard']['collectInsertSizeMetrics']['time'],
        params = config['tools']['picard']['collectInsertSizeMetrics']['params']
    benchmark:
        '{sample}.bam_stats/insert_size_metrics.txt.benchmark'
    threads:
        config['tools']['picard']['collectInsertSizeMetrics']['threads']
    shell:
        ('{config[tools][picard][call]} CollectInsertSizeMetrics ' +
        'I={input.bam} ' +
        'O={output.txt} ' +
        'H={output.pdf} ' +
        '{params.params}')

# This rule can be used to obtain statistics about a sequencing file using fastqc
rule fastqc:
    input:
        fastq = '{fastq}.fastq.gz'
    output:
        html = '{fastq}_fastqc.html',
        zip = '{fastq}_fastqc.zip'
    params:
        lsfoutfile = '{fastq}_fastqc.html.lsfout.log',
        lsferrfile = '{fastq}_fastqc.html.lsferr.log',
        scratch = config['tools']['fastqc']['scratch'],
        mem = config['tools']['fastqc']['mem'],
        time = config['tools']['fastqc']['time'],
        params = config['tools']['fastqc']['params']
    benchmark:
        '{fastq}_fastqc.html.benchmark'
    threads:
        config['tools']['fastqc']['threads']
    shell:
        '{config[tools][fastqc][call]} {params.params} {input.fastq}'

# This rule creates a graphical representation of the relatedness of the samples
if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = BASERECALIBRATIONOUT
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = OUTDIR + 'variants/GATK_HC/raw/'
rule snpHeatmap:
    input:
        vcf = HAPLOTYPECALLEROUT + 'combined.vcf',
    output:
        tsv = temp(HAPLOTYPECALLEROUT + 'combined_dist.tsv'),
        pdf = HAPLOTYPECALLEROUT + 'combined_dist.pdf',
    params:
        lsfoutfile = HAPLOTYPECALLEROUT + 'combined_dist.pdf.lsfout.log',
        lsferrfile = HAPLOTYPECALLEROUT + 'combined_dist.pdf.lsferr.log',
        scratch = config['tools']['snpHeatmap']['scratch'],
        mem = config['tools']['snpHeatmap']['mem'],
        time = config['tools']['snpHeatmap']['time'],
        prepare = config['tools']['snpHeatmap']['prepare']
    benchmark:
        HAPLOTYPECALLEROUT + 'combined_dist.pdf.benchmark'
    threads:
        config['tools']['snpHeatmap']['threads']
    shell:
        ('{params.prepare} {input.vcf} {output.tsv}; ' + 
        '{config[tools][snpHeatmap][call]} {output.tsv} {output.pdf}')
