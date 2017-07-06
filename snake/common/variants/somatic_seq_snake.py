






callerNames = {'varscan_somatic':'VarScan2', 'jointSNVMix2_075': 'JointSNVMix2', 'somaticSniper': 'SomaticSniper', 'varDict': 'VarDict', 'mutect1' : 'MuTect'}
def getModifyMethod(wildcards):
    return callerNames[wildcards.caller]


rule jsm2vcf:
    input:
        jsmfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_classify.jsm'
    output:
        vcffile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}.vcf'
    params:
        lsfoutfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_jsm2vcf.lsfout.log',
        lsferrfile = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_jsm2vcf.lsferr.log',
        scratch = config['tools']['somaticseq']['jsm2vcf']['scratch'],
        mem = config['tools']['somaticseq']['jsm2vcf']['mem'],
        time = config['tools']['somaticseq']['jsm2vcf']['time']
    benchmark:
        JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_jsm2vcf.benchmark'
    shell:
        '{config[tools][somaticseq][jsm2vcf][call]} {input.jsmfile} > {output.vcffile}'

ruleorder: modifymutect > modifyVJSD
rule modifymutect:
    input:
        mutectfile = MUTECT1OUT + '{tumor}_vs_{normal}.vcf',
        tbam = MUTECT1IN + '{tumor}.bam',
        nbam = MUTECT1IN + '{normal}.bam'
    output:
        mutectmodified = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.mutect1_modified.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.mutect1_modify.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.mutect1_modify.lsferr.log',
        scratch = config['tools']['somaticseq']['modifymutect']['scratch'],
        mem = config['tools']['somaticseq']['modifymutect']['mem'],
        time = config['tools']['somaticseq']['modifymutect']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.mutect1_modify.benchmark'
    shell:
        "{config[tools][somaticseq][modifymutect][call]} -type snp -infile {input.mutectfile} -outfile {output.mutectmodified} -nbam {input.nbam} -tbam {input.tbam}"

localrules: symbolicLinkVarScan
rule symbolicLinkVarScan:
    input:
        vcf = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.snp.vcf'
    output:
        vcf = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.vcf'
    params:
        dirName = VARSCANSOMATICOUT
    shell:
        'cd {params.dirName}; ln -s {wildcards.tumor}_vs_{wildcards.normal}.snp.vcf {wildcards.tumor}_vs_{wildcards.normal}.vcf'

ruleorder: symbolicLinkVarDict > modifyVarDict > modifyVJSD
rule modifyVarDict:
    input:
        vcf = VARDICTOUT + '{tumor}_vs_{normal}.vcf'
    output:
        modifiedvcf = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.varDict_modified.vcf',
        modifiedindel = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varDict_modified.vcf'
    params:
        method = 'VarDict',
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modify.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modify.lsferr.log',
        output = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modified.vcf',
        scratch = config['tools']['somaticseq']['modifyvjsd']['scratch'],
        mem = config['tools']['somaticseq']['modifyvjsd']['mem'],
        time = config['tools']['somaticseq']['modifyvjsd']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modify.benchmark'
    shell:
        "{config[tools][somaticseq][modifyvjsd][call]} -filter paired -method {params.method} -infile {input.vcf} -outfile {params.output}"

rule modifyVJSD:
    input:
        vcf = OUTDIR + 'variants/{caller}/raw/{tumor}_vs_{normal}.vcf'
    output:
        modifiedvcf = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.{caller}_modified.vcf'
    params:
        method = getModifyMethod,
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.{caller}_modify.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.{caller}_modify.lsferr.log',
        scratch = config['tools']['somaticseq']['modifyvjsd']['scratch'],
        mem = config['tools']['somaticseq']['modifyvjsd']['mem'],
        time = config['tools']['somaticseq']['modifyvjsd']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.{caller}_modify.benchmark'
    shell:
        "{config[tools][somaticseq][modifyvjsd][call]} -method {params.method} -infile {input.vcf} -outfile {output.modifiedvcf}"

rule modifyVarScanIndel:
    input:
        vcf = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.indel.vcf'
    output:
        modifiedvcf = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varscan_somatic_modified.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varscan_somatic_modify.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varscan_somatic_modify.lsferr.log',
        scratch = config['tools']['somaticseq']['modifyvjsd']['scratch'],
        mem = config['tools']['somaticseq']['modifyvjsd']['mem'],
        time = config['tools']['somaticseq']['modifyvjsd']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varscan_somatic_modify.benchmark'
    shell:
        "{config[tools][somaticseq][modifyvjsd][call]} -method VarScan2 -infile {input.vcf} -outfile {output.modifiedvcf}"

localrules: symbolicLinkVarDict
rule symbolicLinkVarDict:
    input:
        vcf = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.varDict_modified.vcf',
    output:
        vcf = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modified.vcf',
    params:
        dirName = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep'
    shell:
        'cd {params.dirName}; ln -s snp.{wildcards.tumor}_vs_{wildcards.normal}.varDict_modified.vcf {wildcards.tumor}_vs_{wildcards.normal}.varDict_modified.vcf'

rule sortVcfs:
    input:
        modifiedvcf = '{sample}.vcf',
        dict = config['resources'][ORGANISM]['referenceDict']
    output:
        sortedvcf = '{sample}_sorted.vcf'
    params:
        lsfoutfile = '{sample}_sorted.lsfout.log',
        lsferrfile = '{sample}_sorted.lsferr.log',
        scratch = config['tools']['picard']['sortVCF']['scratch'],
        mem = config['tools']['picard']['sortVCF']['mem'],
        time = config['tools']['picard']['sortVCF']['time']
    benchmark:
        '{sample}_sorted.benchmark'
    shell:
        ('{config[tools][picard][call]} SortVcf ' +
        'I={input.modifiedvcf} O={output.sortedvcf} SEQUENCE_DICTIONARY={input.dict}')

rule combineVariantsSNV:
    input:
        mutect1 = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.mutect1_modified_sorted.vcf',
        varDict = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varDict_modified_sorted.vcf',
        somaticSniper = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.somaticSniper_modified_sorted.vcf',
        varScan = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.varscan_somatic_modified_sorted.vcf',
        jsvm2 = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{tumor}_vs_{normal}.jointSNVMix2_075_modified_sorted.vcf',
        hgref = config['resources'][ORGANISM]['reference']
    output:
        combinedVCF = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.combined.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}_combinegatk.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}_combinegatk.lsferr.log',
        scratch = config['tools']['GATK']['combineVariants']['scratch'],
        mem = config['tools']['GATK']['combineVariants']['mem'],
        time = config['tools']['GATK']['combineVariants']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}_combinegatk.benchmark'
    shell:
        ('{config[tools][GATK][call]} -T CombineVariants ' +
        '-R {input.hgref} ' +
        '-nt 1 ' +
        '--genotypemergeoption UNSORTED ' +
        '--variant {input.varScan} ' +
        '--variant {input.somaticSniper} ' +
        '--variant {input.mutect1} ' +
        '--variant {input.jsvm2} ' +
        '--variant {input.varDict} ' +
        '--out {output.combinedVCF} ' +
        '-U ALLOW_SEQ_DICT_INCOMPATIBILITY')

rule combineVariantsIndel:
    input:
        varDict = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varDict_modified_sorted.vcf',
        varScan = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.varscan_somatic_modified_sorted.vcf',
        hgref = config['resources'][ORGANISM]['reference']
    output:
        combinedVCF = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.combined.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}_combinegatk.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}_combinegatk.lsferr.log',
        scratch = config['tools']['GATK']['combineVariants']['scratch'],
        mem = config['tools']['GATK']['combineVariants']['mem'],
        time = config['tools']['GATK']['combineVariants']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}_combinegatk.benchmark'
    shell:
        ('{config[tools][GATK][call]} -T CombineVariants ' +
        ' -R {input.hgref} ' +
        '-nt 1 ' +
        '--genotypemergeoption UNSORTED ' +
        '--variant {input.varScan} ' +
        '--variant {input.varDict} ' +
        '--out {output.combinedVCF} ' +
        '-U ALLOW_SEQ_DICT_INCOMPATIBILITY')

rule scoreSnvs:
    input:
        db = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.combined.dbSNP.cosmic.snpEff.vcf'
    output:
        somaticSNP = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.somaticSNP.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.somaticSNP.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.somaticSNP.lsferr.log',
        scratch = config['tools']['somaticseq']['scoresomaticvariants']['scratch'],
        mem = config['tools']['somaticseq']['scoresomaticvariants']['mem'],
        time = config['tools']['somaticseq']['scoresomaticvariants']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.somaticSNP.benchmark'
    shell:
        '{config[tools][somaticseq][scoresomaticvariants][call]} -tools CGA VarScan2 JointSNVMix2 SomaticSniper VarDict -infile {input.db} -mincaller 1 -outfile {output.somaticSNP}'

rule vcf2tsvSnv:
    input:
        somaticSNP = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.somaticSNP.vcf',
        mutect1Vcf = MUTECT1OUT + '{tumor}_vs_{normal}_sorted.vcf',
        varDictVcf = VARDICTOUT + '{tumor}_vs_{normal}_sorted.vcf',
        somaticSniperVcf = SOMATICSNIPEROUT + '{tumor}_vs_{normal}_sorted.vcf',
        varScanVcf = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.vcf',
        jsvm2Vcf = JOINTSNVMIX2_075_OUT + '{tumor}_vs_{normal}_sorted.vcf',
        hgref = config['resources'][ORGANISM]['reference'],
        #hgref = hg19,
        tbam = JOINTSNVMIX2_075_IN + '{tumor}.bam',
        nbam = JOINTSNVMIX2_075_IN + '{normal}.bam'
    output:
        ensemblessnv = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.Ensemble.sSNV.tsv'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.Ensemble.sSNV.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.Ensemble.sSNV.lsferr.log',
        scratch = config['tools']['somaticseq']['vcf2tsv']['scratch'],
        mem = config['tools']['somaticseq']['vcf2tsv']['mem'],
        time = config['tools']['somaticseq']['vcf2tsv']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/snp.{tumor}_vs_{normal}.Ensemble.sSNV.benchmark'
    shell:
        ('{config[tools][somaticseq][vcf2tsv][call]} -ref {input.hgref} -myvcf {input.somaticSNP} ' +
        '-varscan {input.varScanVcf} ' +
        '-jsm {input.jsvm2Vcf} ' +
        '-sniper {input.somaticSniperVcf} ' +
        '-vardict {input.varDictVcf} ' +
        '-mutect {input.mutect1Vcf} ' +
        '-tbam {input.tbam} -nbam {input.nbam} -dedup -outfile {output.ensemblessnv}')

rule vcf2tsvIndel:
    input:
        somaticSNP = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.somaticSNP.vcf',
        varDictVcf = VARDICTOUT + '{tumor}_vs_{normal}_sorted.vcf',
        varScanVcf = VARSCANSOMATICOUT + '{tumor}_vs_{normal}.vcf',
        hgref = config['resources'][ORGANISM]['reference'],
        #hgref = hg19,
        tbam = VARDICTIN + '{tumor}.bam',
        nbam = VARDICTIN + '{normal}.bam'
    output:
        ensemblessnv = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.Ensemble.sSNV.tsv'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.Ensemble.sSNV.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.Ensemble.sSNV.lsferr.log',
        scratch = config['tools']['somaticseq']['vcf2tsv']['scratch'],
        mem = config['tools']['somaticseq']['vcf2tsv']['mem'],
        time = config['tools']['somaticseq']['vcf2tsv']['mem']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/indel.{tumor}_vs_{normal}.Ensemble.sSNV.benchmark'
    shell:
        ('{config[tools][somaticseq][vcf2tsv][call]} -ref {input.hgref} -myvcf {input.somaticSNP} ' +
        '-varscan {input.varScanVcf} ' +
        '-vardict {input.varDictVcf} ' +
        '-tbam {input.tbam} -nbam {input.nbam} -dedup -outfile {output.ensemblessnv}')

def getAdaClassifier(wildcards):
    if wildcards.type == "indel":
        return config['resources']['general']['adaclassifier_indel']
    elif wildcards.type == "snp":
        return config['resources']['general']['adaclassifier_snp']

rule adaclassify:
    input:
        ensemblessnv = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.tsv',
        classifier = getAdaClassifier
    output:
        trainedSNV = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.trained.tsv'
    params:
        lsfoutfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.trained.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.trained.lsferr.log',
        scratch = config['tools']['somaticseq']['adaclassify']['scratch'],
        mem = config['tools']['somaticseq']['adaclassify']['mem'],
        time = config['tools']['somaticseq']['adaclassify']['time']
    benchmark:
        SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.trained.benchmark'
    shell:
        'R --no-save --args {input.classifier} {input.ensemblessnv} {output.trainedSNV} < {config[tools][somaticseq][adaclassify][call]}'

rule tsv2vcf:
    input:
        trainedSNP = SOMATICSEQOUT + '{tumor}_vs_{normal}_prep/{type}.{tumor}_vs_{normal}.Ensemble.sSNV.trained.tsv'
    output:
        vcf = SOMATICSEQOUT + '{type}.{tumor}_vs_{normal}.somaticseq.vcf'
    params:
        lsfoutfile = SOMATICSEQOUT + '{type}.{tumor}_vs_{normal}.lsfout.log',
        lsferrfile = SOMATICSEQOUT + '{type}.{tumor}_vs_{normal}.lsferr.log',
        scratch = config['tools']['somaticseq']['tsv2vcf']['scratch'],
        mem = config['tools']['somaticseq']['tsv2vcf']['mem'],
        time = config['tools']['somaticseq']['tsv2vcf']['time']
    benchmark:
        SOMATICSEQOUT + '{type}.{tumor}_vs_{normal}.somaticseq.benchmark'
    shell:
        '{config[tools][somaticseq][tsv2vcf][call]} -tsv {input.trainedSNP} -vcf {output.vcf} -tools CGA VarScan2 JointSNVMix2 SomaticSniper VarDict -all'
