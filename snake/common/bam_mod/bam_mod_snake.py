import ntpath

# This rule creates an indes of a BAM file
rule samtools_create_index:
    input:
        bam = '{sample}.bam',
    output:
        idx = '{sample}.bai',
    params:
        lsfoutfile = '{sample}.bai.lsfout.log',
        lsferrfile = '{sample}.bai.lsferr.log',
        scratch = config['tools']['samtools']['index']['scratch'],
        mem = config['tools']['samtools']['index']['mem'],
        time = config['tools']['samtools']['index']['time']
    threads:
        config['tools']['samtools']['index']['threads']
    benchmark:
        '{sample}.bai.benchmark'
    shell:
        '{config[tools][samtools][call]} index {input.bam} && mv {input.bam}.bai {output.idx}'

# This rule create an symbolic link to an existing index
rule linkIndex:
    input:
        bai = '{sample}.bai',
    output:
        bai = '{sample}.bam.bai',
    params:
        lsfoutfile = '{sample}.bam.bai.lsfout.log',
        lsferrfile = '{sample}.bam.bai.lsferr.log',
        scratch = '1000', 
        mem = '1000', 
        time = '1' 
    threads:
        1
    benchmark:
        '{sample}.bam.bai.benchmark'
    shell:
        'dirName=$(dirname "{input.bai}"); inBai=$(basename "{input.bai}"); outBai=$(basename "{output.bai}"); cd "$dirName"; ln -s "$inBai" "$outBai"'

# This rule sorts a BAM file and fixes mate pair information if necessary
if not 'FIXMATEANDSORTIN' in globals():
    FIXMATEANDSORTIN = BWAOUT
if not 'FIXMATEANDSORTOUT' in globals():
    FIXMATEANDSORTOUT = OUTDIR + 'fix_sorted/'
rule picards_fix_mate_pair_and_sort:
    input:
        bam=FIXMATEANDSORTIN + '{sample}.bam'
    output:
        bam=temp(FIXMATEANDSORTOUT + '{sample}.bam')
    params:
        lsfoutfile = FIXMATEANDSORTOUT + '{sample}.bam.lsfout.log',
        lsferrfile = FIXMATEANDSORTOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['picard']['fixMateInformation']['scratch'],
        mem = config['tools']['picard']['fixMateInformation']['mem'],
        time = config['tools']['picard']['fixMateInformation']['time'],
        sortOrder = config['tools']['picard']['fixMateInformation']['sortOrder'],
        assume_sorted = config['tools']['picard']['fixMateInformation']['assume_sorted'],
        params = config['tools']['picard']['fixMateInformation']['params']
    benchmark:
        FIXMATEANDSORTOUT + '{sample}.bam.benchmark'
    threads:
        config['tools']['picard']['fixMateInformation']['threads']
    log:
        FIXMATEANDSORTOUT + '{sample}.bam.log'
    shell:
        ('{config[tools][picard][call]} FixMateInformation ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam} ' +
        'SORT_ORDER={params.sortOrder} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        '2> {log}')

# This functiom creates a list of BAM files created by the read mapper
def getAlignerBams():
    out = []
    if config['tools']['picard']['mergeBams']['useOrphans'] != "Y":
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
    else:
        for f in PAIREDFASTQFILES:
            out.append(os.path.dirname(f) +'/ORPHAN/' + ntpath.basename(f) + '.bam')
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
    return out

# This function is a helper function to get the BAM file names which are then merged.
def getBamsToMerge(wildcards):
    out = []
    allBams = getAlignerBams()
    for bam in allBams:
        if wildcards.sample == bam.split("/")[0]: 
            out.append(MERGEBAMSIN + bam)
    if not out:
        #print(wildcards)
        return ['ERROR']
    return out

# This function prepends 'INPUT=' in front of every BAM file that is to be merged.
def prependBamsToMerge(wildcards):
    bamsToMerge = getBamsToMerge(wildcards)
    return ''.join(['INPUT='+bam+' ' for bam in bamsToMerge])

# This rule merges different BAM files.
if not 'MERGEBAMSIN' in globals():
    MERGEBAMSIN = FIXMATEANDSORTOUT
if not 'MERGEBAMSOUT' in globals():
    MERGEBAMSOUT = OUTDIR + 'merged/'
rule picard_merge_bams:
    input:
        bams = getBamsToMerge
    output:
        bam = temp(MERGEBAMSOUT + '{sample}.bam')
    params:
        lsfoutfile = MERGEBAMSOUT + '{sample}.bam.lsfout.log',
        lsferrfile = MERGEBAMSOUT + '{sample}.bam.lsferr.log',
        assume_sorted = config['tools']['picard']['mergeBams']['assume_sorted'],
        scratch = config['tools']['picard']['mergeBams']['scratch'],
        mem = config['tools']['picard']['mergeBams']['mem'],
        time = config['tools']['picard']['mergeBams']['time'],
        params = config['tools']['picard']['mergeBams']['params'],
        input = prependBamsToMerge
    threads:
        config['tools']['picard']['mergeBams']['threads']
    benchmark:
        MERGEBAMSOUT + '{sample}.bam.benchmark'
    log:
        log = MERGEBAMSOUT + '{sample}.bam.log',
        metrics = MERGEBAMSOUT + '{sample}.bam.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MergeSamFiles ' +
        '{params.input} ' +
        'OUTPUT={output.bam} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        '{params.params} ' +
        '2> {log.log}')

# This rule removes all alignments which are marked as secondary alignment
if not 'NOSECONDARYALNIN' in globals():
    NOSECONDARYALNIN = MERGEBAMSOUT
if not 'NOSECONDARYALNOUT' in globals():
    NOSECONDARYALNOUT = OUTDIR + 'noSecondaryAln/'
rule samtools_remove_secondary_alignments:
    input:
        bam = NOSECONDARYALNIN+ '{sample}.bam',
    output:
        bam=temp(NOSECONDARYALNOUT + '{sample}.bam'),
    params:
        lsfoutfile = NOSECONDARYALNOUT + '{sample}.bam.lsfout.log',
        lsferrfile = NOSECONDARYALNOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['samtools']['rmSecondary']['scratch'],
        mem = config['tools']['samtools']['rmSecondary']['mem'],
        time = config['tools']['samtools']['rmSecondary']['time']
    benchmark:
        NOSECONDARYALNOUT + '{sample}.bam.benchmark'
    threads:
        config['tools']['samtools']['rmSecondary']['threads']
    shell:
        '{config[tools][samtools][call]} view -bh -F 256 {input.bam} > {output.bam}'

# This rule markes PCR duplicates
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = NOSECONDARYALNOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markedDuplicates/'
rule picards_mark_PCR_duplicates:
    input:
        bam=MARKPCRDUBLICATESIN + '{sample}.bam',
    output:
        bam=temp(MARKPCRDUBLICATESOUT + '{sample}.bam'),
    params:
        lsfoutfile = MARKPCRDUBLICATESOUT + '{sample}.bam.lsfout.log',
        lsferrfile = MARKPCRDUBLICATESOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['picard']['markduplicates']['scratch'],
        mem = config['tools']['picard']['markduplicates']['mem'],
        time = config['tools']['picard']['markduplicates']['time'],
        params = config['tools']['picard']['markduplicates']['params']
    threads:
        config['tools']['picard']['markduplicates']['threads']
    benchmark:
        MARKPCRDUBLICATESOUT + '{sample}.bam.benchmark'
    log:
        log = MARKPCRDUBLICATESOUT + '{sample}.bam.log',
        metrics = MARKPCRDUBLICATESOUT + '{sample}.bam.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MarkDuplicates ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam}  ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        'METRICS_FILE={log.metrics} ' +
        '2> {log.log}')

# This rule removed the previously marked PCR duplicates
if not 'REMOVEPCRDUBLICATESIN' in globals():
    REMOVEPCRDUBLICATESIN = MARKPCRDUBLICATESOUT
if not 'REMOVEPCRDUBLICATESOUT' in globals():
    REMOVEPCRDUBLICATESOUT = OUTDIR + 'removedPcrDuplicates/'
rule samtools_remove_PCR_duplicates:
    input:
        bam=REMOVEPCRDUBLICATESIN + '{sample}.bam',
    output:
        bam=REMOVEPCRDUBLICATESOUT + '{sample}.bam',
    params:
        lsfoutfile = REMOVEPCRDUBLICATESOUT + '{sample}.bam.lsfout.log',
        lsferrfile = REMOVEPCRDUBLICATESOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['samtools']['rmDuplicates']['scratch'],
        mem = config['tools']['samtools']['rmDuplicates']['mem'],
        time = config['tools']['samtools']['rmDuplicates']['time']
    threads:
        config['tools']['samtools']['rmDuplicates']['threads']
    benchmark:
        REMOVEPCRDUBLICATESOUT + '{sample}.bam.benchmark'
    shell:
        '{config[tools][samtools][call]} view -bh -F 0x400 {input.bam} > {output.bam}'

# Rule to replace a certain mapping quality with a specified one
# This is a GATK tool
if not 'REASSIGNONEMAPPINGQUALIN' in globals():
    REASSIGNONEMAPPINGQUALIN = OUTDIR + '.reassing_one_mapping_quality_in'
if not 'REASSIGNONEMAPPINGQUALOUT' in globals():
    REASSIGNONEMAPPINGQUALOUT = OUTDIR + '.reassing_one_mapping_quality_out'
rule gatk_reassign_one_mapping_quality_filter:
    input:
        bam=REASSIGNONEMAPPINGQUALIN + '{sample}.bam',
        bamIdx=REASSIGNONEMAPPINGQUALIN + '{sample}.bam.index',
    output:
        bam=REASSIGNONEMAPPINGQUALOUT + '{sample}.bam',
    params:
        scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        time = config['tools']['GATK']['baseRecalibrator']['time'],
        lsfoutfile = REASSIGNONEMAPPINGQUALOUT + '{sample}.bam.lsfout.log',
        lsferrfile = REASSIGNONEMAPPINGQUALOUT + '{sample}.bam.lsferr.log'
    benchmark:
        REASSIGNONEMAPPINGQUALOUT + '{sample}.bam.benchmark'
    threads:
        config['tools']['GATK']['baseRecalibrator']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T PrintReads ' +
        '-R {config[tools][reference]} ' +
        '-I {input.bam} ' +
        '-o {output.bam} ' +
        '-rf ReassignOneMappingQuality ' +
        '-RMQF {config[tools][GATK][reassignOneMappingQualityFilter][oldQual]} ' +
        '-RMQT {config[tools][GATK][reassignOneMappingQualityFilter][newQual]}')

def getSamplesFromExperimentId(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    expMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == 'Y':
                    if exp not in expMap.keys():
                        expMap[exp] = []
                    expMap[exp].append(sample)
                elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
                    if sample in expMap.keys():
                        raise ValueError(sample = " is not uniq in the sample mapping file.")
                    expMap[sample] = []
                    expMap[sample].append(sample)
                else:
                    return "Unknown parameter " + config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] + " to specify whether all bams of one experiment should be realiged together."

    if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "Y":
        if wildcards.experiment not in expMap.keys():
            #print("Experiment " + wildcards.experiment + " is unknown!")
            return "UnknownExperiment"
    elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
        if wildcards.experiment not in expMap.keys():
            #raise ValueError(wildcards.experiment + " is not a valid sample name!")
            return "UnknownSample"
    return expMap[wildcards.experiment]

def getBamsFromExperimentId(wildcards):
    #print("wildcards: ", wildcards)
    #print("getSamplesFromExperimentId: ", getSamplesFromExperimentId(wildcards))
    return expand('{sample}.bam', sample = getSamplesFromExperimentId(wildcards))

def getBaisFromExperimentId(wildcards):
    return expand('{sample}.bai', sample = getSamplesFromExperimentId(wildcards))

def getBamsToRealingFromExperimentId(wildcards):
    #print("wildcards: ", wildcards)
    #print("getBamsFromExperimentId: ", getBamsFromExperimentId(wildcards))
    return expand(REALIGNINDELSIN + '{bam}', bam = getBamsFromExperimentId(wildcards))

def getBaisToRealingFromExperimentId(wildcards):
    return expand(REALIGNINDELSIN + '{bai}', bai = getBaisFromExperimentId(wildcards))

def prependBamsToRealign(wildcards):
    bamsToRealign = getBamsToRealingFromExperimentId(wildcards)
    return ''.join(['-I '+bam+' ' for bam in bamsToRealign])

def getDataBasisForRealign():
    out = []
    out.append(config['resources'][ORGANISM]['reference']) # this is a dummy such that something is retured
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['Mills_indels'])
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['1000G_indels'])
    return out

def prependDataBasisForRealignTargetCreator():
    out = ""
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['Mills_indels'], Error):
            print("You have not specified config[resources][ORGANISM][Mills_indels]")
        out += " --known " + config['resources'][ORGANISM]['Mills_indels']
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['1000G_indels'], Error):
            print("You have not specified config[resources][ORGANISM][1000G_indels]")
        out += " --known " + config['resources'][ORGANISM]['1000G_indels']
    return out

# Rule to create a file containing the regions to perfrom the indel realignment on
if not 'REALIGNINDELSIN' in globals():
    REALIGNINDELSIN = REMOVEPCRDUBLICATESOUT
if not 'REALIGNINDELSOUT' in globals():
    REALIGNINDELSOUT = OUTDIR + 'realignedIndels/'
rule gatk_realign_target_creation:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId,
        reference = config['resources'][ORGANISM]['reference'],
        databasis = getDataBasisForRealign()
    output:
        intervals = temp(REALIGNINDELSOUT + '{experiment}.intervals'),
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}.intervals.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}.intervals.lsferr.log',
        scratch = config['tools']['GATK']['realign']['targetCreator']['scratch'],
        mem = config['tools']['GATK']['realign']['targetCreator']['mem'],
        time = config['tools']['GATK']['realign']['targetCreator']['time'],
        params = config['tools']['GATK']['realign']['targetCreator']['params'],
        input = prependBamsToRealign,
        known = prependDataBasisForRealignTargetCreator(),
    benchmark:
        REALIGNINDELSOUT + '{experiment}.intervals.benchmark'
    threads:
        config['tools']['GATK']['realign']['targetCreator']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T RealignerTargetCreator ' +
        '-R {input.reference} ' +
        '{params.input} ' +
        '{params.known} ' +
        '-o {output.intervals} ' +
        '-nt {threads} ' +
        '{params.params}')

localrules: createRealingIndelsInOutMapping
rule createRealingIndelsInOutMapping:
    output:
        exp = REALIGNINDELSOUT + '{experiment}.map'
    run:
        import sys
        try:
            outFile = open(output.exp, 'x')
            bams = getBamsFromExperimentId(wildcards)
            for bam in bams:
                outFile.write(bam + "\t" + REALIGNINDELSOUT + "ORIGINAL_" + bam + "\n")
        except IOError:
            print("Could not open file: ", output.exp)

def prependDataBasisForTargetRealigner():
    out = ""
    if config['tools']['GATK']['realign']['Mills_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['Mills_indels'], Error):
            print("You have not specified config[resources][ORGANISM][Mills_indels]")
        out += " -known " + config['resources'][ORGANISM]['Mills_indels']
    if config['tools']['GATK']['realign']['1000G_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['1000G_indels'], Error):
            print("You have not specified config[resources][ORGANISM][1000G_indels]")
        out += " -known " + config['resources'][ORGANISM]['1000G_indels']
    return out
# Rule to perform the indel realignment
# This is a GATK tool
#ruleorder: realignIndels > createIndex
rule gatk_realign_indels:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId,
        map = REALIGNINDELSOUT + '{experiment}.map',
        intervals = REALIGNINDELSOUT + '{experiment}.intervals',
        reference = config['resources'][ORGANISM]['reference'],
        databasis = getDataBasisForRealign()
    output:
        txt = REALIGNINDELSOUT + '{experiment}.realigned.txt',
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}.realigned.txt.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}.realigned.txt.lsferr.log',
        scratch = config['tools']['GATK']['realign']['realignIndels']['scratch'],
        mem = config['tools']['GATK']['realign']['realignIndels']['mem'],
        time = config['tools']['GATK']['realign']['realignIndels']['time'],
        params = config['tools']['GATK']['realign']['realignIndels']['params'],
        input = prependBamsToRealign,
        known = prependDataBasisForTargetRealigner(),
    benchmark:
        REALIGNINDELSOUT + '{experiment}.realigned.txt.benchmark'
    threads:
        config['tools']['GATK']['realign']['realignIndels']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T IndelRealigner ' +
        '-R {input.reference} ' +
        '{params.input} ' +
        '{params.known} ' +
        '--nWayOut {input.map} ' +
        '-targetIntervals {input.intervals} ' +
        '{params.params} ' +
        '&& touch {output.txt}')

def getExperimentIdFromBam(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "Y":
                    sampleMap[sample] = exp
                elif config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] == "N":
                    sampleMap[sample] = sample
                else:
                    return "Unknown parameter " + config['tools']['GATK']['realign']['realignFilesFromExperimentTogether'] + " to specify whether all bams of one experiment should be realiged together."
    
    if wildcards.sample not in sampleMap.keys():
        #raise ValueError(wildcards.sample + " is not a sample ID")
        return "UnknownSample"
    return sampleMap[wildcards.sample]

def getRealignedExperiment(wildcards):
    return expand(REALIGNINDELSOUT + '{experiment}.realigned.txt', experiment = getExperimentIdFromBam(wildcards))


# This is a "dummy" rule
localrules: getRealignedBam 
rule getRealignedBam:
    input:
        txt = getRealignedExperiment
    output:
        bam = REALIGNINDELSOUT + '{sample}.bam'
    params:
        dirName = REALIGNINDELSOUT,
        originalBam = REALIGNINDELSOUT + 'ORIGINAL_{sample}.bam'
    shell:
        'cd {params.dirName}; ln -s ORIGINAL_{wildcards.sample}.bam {wildcards.sample}.bam'

def getDataBasisForBaseRecalibration():
    out = []
    out.append(config['resources'][ORGANISM]['reference']) # this is a dummy such that something is retured
    if config['tools']['GATK']['baseRecalibrator']['Mills_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['Mills_indels'])
    if config['tools']['GATK']['baseRecalibrator']['1000G_indels'] == "Y":
        out.append(config['resources'][ORGANISM]['1000G_indels'])
    if config['tools']['GATK']['baseRecalibrator']['dbSNP'] == "Y":
        out.append(config['resources'][ORGANISM]['dbSNP'])
    return out

def prependDataBasisForBaseRecalibration():
    out = ""
    if config['tools']['GATK']['baseRecalibrator']['Mills_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['Mills_indels'], Error):
            print("You have not specified config[resources][ORGANISM][Mills_indels]")
        out += " -knownSites " + config['resources'][ORGANISM]['Mills_indels']
    if config['tools']['GATK']['baseRecalibrator']['1000G_indels'] == "Y":
        if isinstance(config['resources'][ORGANISM]['1000G_indels'], Error):
            print("You have not specified config[resources][ORGANISM][1000G_indels]")
        out += " -knownSites " + config['resources'][ORGANISM]['1000G_indels']
    if config['tools']['GATK']['baseRecalibrator']['dbSNP'] == "Y":
        if isinstance(config['resources'][ORGANISM]['dbSNP'], Error):
            print("You have not specified config[resources][ORGANISM][dbSNP]")
        out += " -knownSites " + config['resources'][ORGANISM]['dbSNP']
    return out
# Rule to create the baserecalibration table used for the baserecalibration
# This is a GATK tool
if not 'BASERECALIBRATIONIN' in globals():
    BASERECALIBRATIONIN = REALIGNINDELSOUT
if not 'BASERECALIBRATIONOUT' in globals():
    BASERECALIBRATIONOUT = OUTDIR + 'recalibratedBases/'
rule gatk_first_pass_create_recalibration_table:
    input:
        bam = BASERECALIBRATIONIN + '{sample}.bam',
        bai = BASERECALIBRATIONIN + '{sample}.bai',
        reference = config['resources'][ORGANISM]['reference'],
        referenceDict = config['resources'][ORGANISM]['referenceDict'],
        regions = config['resources'][ORGANISM]['regions'],
        databasis = getDataBasisForBaseRecalibration()
    output:
        tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
    params:
        lsfoutfile = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table.lsfout.log',
        lsferrfile = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table.lsferr.log',
        scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        time = config['tools']['GATK']['baseRecalibrator']['time'],
        params = config['tools']['GATK']['baseRecalibrator']['params'],
        known = prependDataBasisForBaseRecalibration()
    benchmark:
        BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table.benchmark'
    threads:
        config['tools']['GATK']['baseRecalibrator']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T BaseRecalibrator ' +
        '-L {input.regions} ' +
        '-R {input.reference} ' +
        '-I {input.bam} ' +
        '{params.known} ' + 
        '-nct {threads} ' +
        '{params.params} ' +        
        '-o {output.tab}')

# Rule to create the base-recalibration table used to analyze the effect of the baserecalibration
# This is a GATK tool
#ruleorder: secondPassCreateRecalibrationTable > firstPassCreateRecalibrationTable
rule gatk_second_pass_create_recalibration_table:
    input:
        bam = BASERECALIBRATIONOUT + '{sample}.bam',
        bai = BASERECALIBRATIONOUT + '{sample}.bai',
        tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
        reference = config['resources'][ORGANISM]['reference'],
        referenceDict = config['resources'][ORGANISM]['referenceDict'],
        regions = config['resources'][ORGANISM]['regions'],
        databasis = getDataBasisForBaseRecalibration()
    output:
        tab = BASERECALIBRATIONOUT + '{sample}.secondPass_reca.table',
    params:
        lsfoutfile = BASERECALIBRATIONOUT + '{sample}_secondPass_reca.table.lsfout.log',
        lsferrfile = BASERECALIBRATIONOUT + '{sample}_secondPass_reca.table.lsferr.log',
        scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        time = config['tools']['GATK']['baseRecalibrator']['time'],
        params = config['tools']['GATK']['baseRecalibrator']['params'],
        known = prependDataBasisForBaseRecalibration()
    benchmark:
        BASERECALIBRATIONOUT + '{sample}_secondPass_reca.table.benchmark'
    threads:
        config['tools']['GATK']['baseRecalibrator']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T BaseRecalibrator ' +
        '-L {input.regions} ' +
        '-R {input.reference} ' +
        '-I {input.bam} ' +
        '{params.known} ' +
        '-nct {threads} ' +
        '-BQSR {input.tab} ' +
        '{params.params} ' +                
        '-o {output.tab}')

# Rule to realing the reads around indels
# This is a GATK tool
rule gatk_base_recalibration:
    input:
        bam = BASERECALIBRATIONIN + '{sample}.bam',
        bai = BASERECALIBRATIONIN + '{sample}.bai',
        tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
        reference = config['resources'][ORGANISM]['reference'],
        referenceDict = config['resources'][ORGANISM]['referenceDict']
    output:
        bam = BASERECALIBRATIONOUT + '{sample}.bam',
    params:
        lsfoutfile = BASERECALIBRATIONOUT + '{sample}.bam.lsfout.log',
        lsferrfile = BASERECALIBRATIONOUT + '{sample}.bam.lsferr.log',
        scratch = config['tools']['GATK']['baseRecalibratorPrintReads']['scratch'],
        mem = config['tools']['GATK']['baseRecalibratorPrintReads']['mem'],
        params = config['tools']['GATK']['baseRecalibratorPrintReads']['params'],        
        time = config['tools']['GATK']['baseRecalibratorPrintReads']['time']
    benchmark:
        BASERECALIBRATIONOUT + '{sample}.bam.benchmark'
    threads:
        config['tools']['GATK']['baseRecalibratorPrintReads']['threads']
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T PrintReads ' +
        '-R {input.reference} ' +
        '-I {input.bam} ' +
        '-o {output.bam} ' +
        '{params.params} ' +                        
        '-BQSR {input.tab}')

if not 'MPILEUPIN' in globals():
    MPILEUPIN = BASERECALIBRATIONOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'
rule mpileupBcf:
    input:
        bam = MPILEUPIN + '{sample}.bam',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        bcf = temp(MPILEUPOUT + '{sample}.bcf')
    params:
        lsfoutfile = MPILEUPOUT + '{sample}.bcf.lsfout.log',
        lsferrfile = MPILEUPOUT + '{sample}.bcf.lsferr.log',
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
        params = config['tools']['samtools']['mpileup']['params']
    threads:
        config['tools']['samtools']['mpileup']['threads']
    benchmark:
        MPILEUPOUT + '{sample}.bcf.benchmark'
    log:
        MPILEUPOUT + '{sample}.bcf.log'
    shell:
        ('{config[tools][samtools][call]} mpileup ' +
        '{params.params} ' + 
        '-f {input.reference} ' +
        '-o {output.bcf} ' +
        '-l {input.regions} ' +
        '{input.bam}')

rule mpileupMpileup:
    input:
        bam = MPILEUPIN + '{sample}.bam',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        mpileup = temp(MPILEUPOUT + '{sample}.mpileup')
    params:
        lsfoutfile = MPILEUPOUT + '{sample}.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{sample}.mpileup.lsferr.log',
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
        params = config['tools']['samtools']['mpileup']['params']
    threads:
        config['tools']['samtools']['mpileup']['threads']
    benchmark:
        MPILEUPOUT + '{sample}.mpileup.benchmark'
    log:
        MPILEUPOUT + '{sample}.mpileup.log'
    shell:
        ('{config[tools][samtools][call]} mpileup ' +
        '{params.params} ' + 
        '-f {input.reference} ' + 
        '-o {output.mpileup} ' + 
        '-l {input.regions} ' +
        '{input.bam}')

