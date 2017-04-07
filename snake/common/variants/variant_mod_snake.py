# This rule uses a python script to fill the header of VarScan and freebayes vcf files
if not 'VARSCANUPDATEHEADERIN' in globals():
    VARSCANUPDATEHEADERIN = VARSCANSOMATICOUT
if not 'VARSCANUPDATEHEADEROUT' in globals():
    VARSCANUPDATEHEADEROUT = OUTDIR + 'variants/varscan_somatic/complete_raw/'
rule updateVCFHeader:
    input:
        vcf = VARSCANUPDATEHEADERIN + '{sample}.vcf',
        reference = config['resources'][ORGANISM]['referenceNamesForVcf']
    output:
        vcf = VARSCANUPDATEHEADEROUT + '{sample}.vcf'
    params:
        lsfoutfile = VARSCANUPDATEHEADEROUT + '{sample}.vcf.lsfout.log',
        lsferrfile = VARSCANUPDATEHEADEROUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['updateVCFHeader']['scratch'],
        mem = config['tools']['updateVCFHeader']['mem'],
        time = config['tools']['updateVCFHeader']['time']
    threads:
        config['tools']['updateVCFHeader']['threads']
    benchmark:
        VARSCANUPDATEHEADEROUT + '{sample}.vcf.benchmark'
    shell:
        '{config[tools][updateVCFHeader][call]} {input.vcf} {input.reference} {output.vcf}'
        
# This rule uses bcftools convert to compress a vcf
rule bcftoolsConvert:
    input:
        vcf = '{sample}.vcf'
    output:
        vcf = '{sample}.vcf.gz'
    params:
        lsfoutfile = '{sample}.vcf.gz.lsfout.log',
        lsferrfile = '{sample}.vcf.gz.lsferr.log',
        scratch = config['tools']['bcftools']['scratch'],
        mem = config['tools']['bcftools']['mem'],
        time = config['tools']['bcftools']['time']
    threads:
        config['tools']['bcftools']['threads']
    benchmark:
        '{sample}.vcf.gz.benchmark'
    shell:
        '{config[tools][bcftools][call]} convert -Oz -o {output.vcf} {input.vcf}'
        
# This rule uses bcftools index to index a compressed vcf
rule bcftoolsIndex:
    input:
        vcf = '{sample}.vcf.gz'
    output:
        vcf = '{sample}.vcf.gz.csi'
    params:
        lsfoutfile = '{sample}.vcf.gz.csi.lsfout.log',
        lsferrfile = '{sample}.vcf.gz.csi.lsferr.log',
        scratch = config['tools']['bcftools']['scratch'],
        mem = config['tools']['bcftools']['mem'],
        time = config['tools']['bcftools']['time']
    threads:
        config['tools']['bcftools']['threads']
    benchmark:
        '{sample}.vcf.gz.csi.benchmark'
    shell:
        '{config[tools][bcftools][call]} index {input.vcf}'
        
# This rule uses bcftools index to index a compressed vcf
if not 'VARSCANCOMPLETEIN' in globals():
    VARSCANCOMPLETEIN = VARSCANUPDATEHEADEROUT
if not 'VARSCANCOMPLETEOUT' in globals():
    VARSCANCOMPLETEOUT = OUTDIR + 'variants/varscan_somatic/combined_raw/'
ruleorder: snpSift_dbNSFP_Annotation > snpSift_COSMIC_Annotation > snpSift_clinVar_Annotation > snpSift_dbSNP_Annotation > snpEffAnnotation > updateNormalTumorName > bcftoolsConcat
rule bcftoolsConcat:
    input:
        vcfIndel = VARSCANCOMPLETEIN + '{sample}.indel.vcf.gz',
        vcfSnp = VARSCANCOMPLETEIN + '{sample}.snp.vcf.gz',
        vcfIndelIndex = VARSCANCOMPLETEIN + '{sample}.indel.vcf.gz.csi',
        vcfSnpIndex = VARSCANCOMPLETEIN + '{sample}.snp.vcf.gz.csi'
    output:
        vcf = VARSCANCOMPLETEOUT + '{sample}.vcf'
    params:
        lsfoutfile = VARSCANCOMPLETEOUT + '{sample}.vcf.lsfout.log',
        lsferrfile = VARSCANCOMPLETEOUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['bcftools']['scratch'],
        mem = config['tools']['bcftools']['mem'],
        time = config['tools']['bcftools']['time']
    threads:
        config['tools']['bcftools']['threads']
    benchmark:
        VARSCANCOMPLETEOUT + '{sample}.vcf.benchmark'
    shell:
        '{config[tools][bcftools][call]} concat {input.vcfSnp} {input.vcfIndel} -a -o {output.vcf}'

# This rule annotates a vcf file using snpEff
ruleorder: snpSift_dbNSFP_Annotation > snpSift_COSMIC_Annotation > snpSift_clinVar_Annotation > snpSift_dbSNP_Annotation > snpEffAnnotation > updateNormalTumorName > mutect1
ruleorder: snpSift_dbNSFP_Annotation > snpSift_COSMIC_Annotation > snpSift_clinVar_Annotation > snpSift_dbSNP_Annotation > snpEffAnnotation > updateNormalTumorName > strelka
# NOTE: all anotation calls could theoretically be combined into one pipe command; however, this could lead to problems when not all databases available!
rule snpEffAnnotation:
    input:
        vcf = '{sample}.vcf',
        snpEffDB  = {config['resources'][ORGANISM]['pathSnpEffDB']}
    output:
        vcf = '{sample}.snpEff.vcf',
        stats = '{sample}.snpEff.stats.html'
    params:
        lsfoutfile = '{sample}.snpEff.vcf.lsfout.log',
        lsferrfile = '{sample}.snpEff.vcf.lsferr.log',
        scratch = config['tools']['snpEff']['scratch'],
        mem = config['tools']['snpEff']['mem'],
        time = config['tools']['snpEff']['time'],
        dbName = config['tools']['snpEff']['dbName']
    threads:
        config['tools']['snpEff']['threads']
    benchmark:
        '{sample}.snpEff.vcf.benchmark'
    shell: 
        '{config[tools][snpEff][call]} ann {params.dbName} -dataDir {input.snpEffDB} -nodownload -s {output.stats} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the dbSNP database
rule snpSift_dbSNP_Annotation:
    input:
        vcf = '{sample}.vcf',
        dbSnpDB  = {config['resources'][ORGANISM]['dbSNP']}
    output:
        vcf = '{sample}.dbSNP.vcf'
    params:
        lsfoutfile = '{sample}.dbSNP.vcf.lsfout.log',
        lsferrfile = '{sample}.dbSNP.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        config['tools']['snpSift']['threads']
    benchmark:
        '{sample}.dbSNP.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.dbSnpDB} {input.vcf} > {output.vcf}'

# This rule annotates a vcf file using snpSift and the clinVar database
rule snpSift_clinVar_Annotation:
    input:
        vcf = '{sample}.vcf',
        clinVarDB  = {config['resources'][ORGANISM]['clinvar']}
    output:
        vcf = '{sample}.clinVar.vcf'
    params:
        lsfoutfile = '{sample}.clinVar.vcf.lsfout.log',
        lsferrfile = '{sample}.clinVar.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        config['tools']['snpSift']['threads']
    benchmark:
        '{sample}.clinVar.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.clinVarDB} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the cosmic database
rule snpSift_COSMIC_Annotation:
    input:
        vcf = '{sample}.vcf',
        cosmicDB  = config['resources'][ORGANISM]['cosmic']
    output:
        vcf = '{sample}.cosmic.vcf'
    params:
        lsfoutfile = '{sample}.cosmic.vcf.lsfout.log',
        lsferrfile = '{sample}.cosmic.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        config['tools']['snpSift']['threads']
    benchmark:
        '{sample}.cosmic.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.cosmicDB} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the dbNSFP database (functional annotation)
rule snpSift_dbNSFP_Annotation:
    input:
        vcf = '{sample}.vcf',
        dbNSFPDB  = config['resources'][ORGANISM]['dbnsfp']
    output:
        vcf = '{sample}.annotated.vcf'
    params:
        lsfoutfile = '{sample}.annotated.vcf.lsfout.log',
        lsferrfile = '{sample}.annotated.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        config['tools']['snpSift']['threads']
    benchmark:
        '{sample}.annotated.vcf.benchmark'
    shell:
        'dbPath=$(readlink {input.dbNSFPDB}) ; {config[tools][snpSift][call]} dbnsfp -db $dbPath {input.vcf} > {output.vcf}'
        
# This rule filters all "REJECT" lines from the mutect1 result
if not 'MUTECT1FILTERIN' in globals():
    MUTECT1FILTERIN = MUTECT1OUT
if not 'MUTECT1FILTEROUT' in globals():
    MUTECT1FILTEROUT = OUTDIR + 'variants/mutect1/filtered/'
ruleorder: filterMutect1Reject > mutect1
rule filterMutect1Reject:
    input:
        vcf = MUTECT1FILTERIN + '{sample}.annotated.vcf',
        txt = MUTECT1FILTERIN + '{sample}.txt'
    output:
        vcf = MUTECT1FILTEROUT + '{sample}.annotated.pass.vcf',
        txt = MUTECT1FILTEROUT + '{sample}.pass.txt'
    params:
        lsfoutfile = MUTECT1FILTEROUT + '{sample}.annotated.pass.vcf.lsfout.log',
        lsferrfile = MUTECT1FILTEROUT + '{sample}.annotated.pass.vcf.lsferr.log',
        scratch = config['tools']['simpleMutect1Filter']['scratch'],
        mem = config['tools']['simpleMutect1Filter']['mem'],
        time = config['tools']['simpleMutect1Filter']['time']
    threads:
        config['tools']['simpleMutect1Filter']['threads']
    benchmark:
        MUTECT1FILTEROUT + '{sample}.annotated.pass.vcf.benchmark'
    shell:
        '{config[tools][simpleMutect1Filter][call]} {input.vcf} {input.txt} {output.vcf} {output.txt}'
        
# This rule filters rejected lines from the strelka result
if not 'STRELKAFILTERIN' in globals():
    STRELKAFILTERIN = STRELKAOUT
if not 'STRELKAFILTEROUT' in globals():
    STRELKAFILTEROUT = OUTDIR + 'variants/strelka/filtered/'
ruleorder: filterStrelka > strelka
rule filterStrelka:
    input:
        vcf = STRELKAFILTERIN + '{sample}.annotated.vcf'
    output:
        vcf = STRELKAFILTEROUT + '{sample}.annotated.pass.vcf'
    params:
        lsfoutfile = STRELKAFILTEROUT + '{sample}.annotated.pass.vcf.lsfout.log',
        lsferrfile = STRELKAFILTEROUT + '{sample}.annotated.pass.vcf.lsferr.log',
        scratch = config['tools']['strelkaFilter']['scratch'],
        mem = config['tools']['strelkaFilter']['mem'],
        time = config['tools']['strelkaFilter']['time']
    threads:
        config['tools']['strelkaFilter']['threads']
    benchmark:
        STRELKAFILTEROUT + '{sample}.annotated.pass.vcf.benchmark'
    shell:
        '{config[tools][strelkaFilter][call]} {input.vcf} {output.vcf}'
        
# This rule filters annotated vcf files produced by VarScan2 somatic
if not 'VARSCANSOMATICFILTERIN' in globals():
    VARSCANSOMATICFILTERIN = VARSCANCOMPLETEOUT
if not 'VARSCANSOMATICFILTEROUT' in globals():
    VARSCANSOMATICFILTEROUT = OUTDIR + 'variants/varscan_somatic/filtered/'
rule filterVarScanSomatic:
    input:
        vcf = VARSCANSOMATICFILTERIN + '{sample}.vcf'
    output:
        vcf = VARSCANSOMATICFILTEROUT + '{sample}.vcf'
    params:
        lsfoutfile = VARSCANSOMATICFILTEROUT + '{sample}.vcf.lsfout.log',
        lsferrfile = VARSCANSOMATICFILTEROUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['varscanSomaticFilter']['scratch'],
        mem = config['tools']['varscanSomaticFilter']['mem'],
        time = config['tools']['varscanSomaticFilter']['time'],
        minVarSupport = config['tools']['varscanSomaticFilter']['minVarSupport'],
        pvalue = config['tools']['varscanSomaticFilter']['pvalue'],
        minNucCoverage = config['tools']['varscanSomaticFilter']['minNucCoverage'],
        filterStrands = config['tools']['varscanSomaticFilter']['filterStrands'],
        tumorFreqThreshold = config['tools']['varscanSomaticFilter']['tumorFreqThreshold'],
        filterHomopolymer = config['tools']['varscanSomaticFilter']['filterHomopolymer'],
        filterSilent = config['tools']['varscanSomaticFilter']['filterSilent'],
        lohThreshold = config['tools']['varscanSomaticFilter']['lohThreshold']
    threads:
        config['tools']['varscanSomaticFilter']['threads']
    benchmark:
        VARSCANSOMATICFILTEROUT + '{sample}.vcf.benchmark'
    shell:
        ('{config[tools][varscanSomaticFilter][call]} {input.vcf} {output.vcf} ' +
        '{params.minVarSupport} ' +
        '{params.pvalue} ' +
        '{params.minNucCoverage} ' +
        '{params.filterStrands} ' +
        '{params.tumorFreqThreshold} ' +
        '{params.filterHomopolymer} ' +
        '{params.filterSilent} ' +
        '{params.lohThreshold}')

# this rule fixes the header of somaticSniper, varScan2, freebayes, strelka... 
# NORMAL and TUMOR are replaced with the correct sample names 
# -> this is done immediately after creating the vcfs
rule updateNormalTumorName:
    input:
        vcf = '{tumor}_vs_{normal}.vcf' 
    output:
        vcf = '{tumor}_vs_{normal}.correctNames.vcf'
    params:
        tumor = '{tumor}',
        normal = '{normal}',
        lsfoutfile = '{tumor}_vs_{normal}.correctNames.vcf.lsfout.log',
        lsferrfile = '{tumor}_vs_{normal}.correctNames.vcf.lsferr.log',
        scratch = config['tools'][ 'updateNormalTumorName']['scratch'],
        mem = config['tools']['updateNormalTumorName']['mem'],
        time = config['tools']['updateNormalTumorName']['time']
    threads:
        config['tools']['updateNormalTumorName']['threads']
    benchmark:
        '{tumor}_vs_{normal}.correctNames.vcf.benchmark'
    shell:
        ('tumorName=$(basename {params.tumor}) ; sed \"/^#CHROM/s/TUMOR/$tumorName/\" {input.vcf} > {output.vcf} ; ' +
        'sed -i \"/^#CHROM/s/NORMAL/{params.normal}/\" {output.vcf}')
        
# This rule combines the variant calls of several callers, using GATK VriantCombine
rule gatk_variantCombine:
    input:
        vcf = '{sample}.annotated.vcf'
    output:
        vcf = '{sample}.combined.vcf'
    params:
        lsfoutfile = '{sample}.annotated.vcf.lsfout.log',
        lsferrfile = '{sample}.annotated.vcf.lsferr.log',
        scratch = config['tools']['GATK']['combineVariants']['scratch'],
        mem = config['tools']['GATK']['combineVariants']['mem'],
        time = config['tools']['GATK']['combineVariants']['time']
    threads:
        config['tools']['GATK']['combineVariants']['threads']
    benchmark:
        '{sample}.annotated.vcf.benchmark'
    shell:
        '{config[tools][GATK][call]} -T -db {input.dbNSFPDB} {input.vcf} > {output.vcf}'
