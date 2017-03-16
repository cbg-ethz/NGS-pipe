import os.path

def getSampleNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
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
                    sample = lineSplit[1]
                    output.append(sample)
    return output

def getSingleFastqFiles():
    return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

def getPairedFastqFiles():
    return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
    return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq')]

def getPairedFastqFilesWithoutR():
    return [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    return [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

def getNormalTumorFiles():
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    output = []
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if exp not in sampleMap.keys():
                    sampleMap[exp] = dict()
                if tpoint not in sampleMap[exp].keys():
                    sampleMap[exp][tpoint] = dict()
                if sampleType not in sampleMap[exp][tpoint].keys():
                    sampleMap[exp][tpoint][sampleType] = []
                sampleMap[exp][tpoint][sampleType].append(sample)
    for expKey, expValue in sampleMap.items():
        for tpointKey, tpointValue in expValue.items():
            if 'T' in tpointValue and 'N' in tpointValue:
                for sampleTumor in tpointValue['T']:
                    for sampleNormal in tpointValue['N']:
                        output.append(sampleTumor + '_vs_' + sampleNormal)
    return output


#def replaceResource(line):
#    out=''
#    splitValues=line.split('$$$')
#    for i in range(0, len(splitValues)):
#        if not i%2:
#         
#            out += splitValues[i]
#        else:
#            location = splitValues[i].split(':')
#            if len(location) == 2:
#                try:
#                    out += config[location[0]][location[1]]
#                except KeyError:
#                    return ''
#            if len(location) == 3:
#                try:
#                    if location[1] == "general":
#                        out += config[location[0]][location[1]][location[2]]
#                    else:
#                        out += config[location[0]][ORGANISM][location[2]]
#                except KeyError:
#                    return ''
#    return out
#
#def replaceResourceRecursive(dictionary):
#    for key, value in dictionary.items():
#        if isinstance(value, dict):
#            replaceResourceRecursive(value)
#        else:
#            dictionary[key] = replaceResource(value)
#
#def postProcessConfigMap():
#    # complete all file path in the dictionary
#    global TOOLSDIR
#    global RESOURCEDIR
#    if TOOLSDIR[-1] != '/':
#        TOOLSDIR = TOOLSDIR + '/'
#    if RESOURCEDIR[-1] != '/':
#        RESOURCEDIR += '/'
#    config['dirs']={"tools": TOOLSDIR, "resource": RESOURCEDIR}
#    for organism in config['resources']:
#        for resource in config['resources'][organism]:
#            if resource in config['resources']['projectSpecific']:
#                config['resources'][organism][resource] = config['resources']['projectSpecific'][resource]
#            else:
#                config['resources'][organism][resource] = config['dirs']['resource'] + organism + '/' + config['resources'][organism][resource]
#    for tool, toolSpecifications in config['tools'].items():
#        replaceResourceRecursive(toolSpecifications)
#
#    # GATK specific
#    # if not all databases are present adapt the config map
#    if not os.path.isfile(config['resources'][ORGANISM]['Mills_indels']):
#        print("WARNING: ", config['resources'][ORGANISM]['Mills_indels'], " not present. GATK RealignTargetCreator and GATK BaseRecalibration will not be able to use it.")
#        config['tools']['GATK']['realignTargetCreator']['known1'] = ''
#        config['tools']['GATK']['baseRecalibrator']['known1'] = ''
#    if not os.path.isfile(config['resources'][ORGANISM]['1000G_indels']):
#        print("WARNING: ", config['resources'][ORGANISM]['1000G_indels'], " not present. GATK RealignTargetCreator and GATK BaseRecalibration will not be able to use it.")
#        config['tools']['GATK']['realignTargetCreator']['known2']= ''
#        config['tools']['GATK']['baseRecalibrator']['known2'] = ''
#    if not os.path.isfile(config['resources'][ORGANISM]['dbSNP']):
#        print("WARNING: ", config['resources'][ORGANISM]['dbSNP'], " not present. GATK BaseRecalibration and GATK SNPRecalibrationModel will not be able to use it.")
#        config['tools']['GATK']['baseRecalibrator']['known3'] = ''
#        config['tools']['GATK']['gatkSNPrecalibrateModel']['resource4'] = ''
#    if not os.path.isfile(config['resources'][ORGANISM]['hapmap']):
#        print("WARNING: ", config['resources'][ORGANISM]['hapmap'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
#        config['tools']['GATK']['gatkSNPrecalibrateModel']['resource1'] = ''
#    if not os.path.isfile(config['resources'][ORGANISM]['1000G_omni']):
#        print("WARNING: ", config['resources'][ORGANISM]['1000G_omni'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
#        config['tools']['GATK']['gatkSNPrecalibrateModel']['resource2'] = ''
#    if not os.path.isfile(config['resources'][ORGANISM]['1000G_snp']):
#        print("WARNING: ", config['resources'][ORGANISM]['1000G_snp'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
#        config['tools']['GATK']['gatkSNPrecalibrateModel']['resource3'] = ''
#    if not os.path.isfile(config['resources']['general']['gatkKey']):
#        print("WARNING: GATK will upload a report for each tool used to the Amazon cloud because no GATK key was provided!")
#        config['tools']['GATK']['call'] = config['tools']['GATK']['call'].strip().split(" --gatk_key")[0]

#rule gunzip:
#    input: 
#        '{sample}.gz'
#    output: 
#        '{sample}'
#    params:
#        lsfoutfile = '{sample}.gunzip.lsfout.log',
#        lsferrfile = '{sample}.gunzip.lsferr.log',
#        scratch = config['tools']['gunzip']['scratch'],
#        mem = config['tools']['gunzip']['memory'],
#        time = config['tools']['gunzip']['time']
#    shell:
#        'gunzip {input}'
