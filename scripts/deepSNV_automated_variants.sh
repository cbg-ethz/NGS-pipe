#!/bin/bash -l

#function deepSNV_automated_variants()
#{
	if [ -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" ]; then #
		echo "usage deepSNV_automated_variants <bed_file> <tumor_bam> <normal_bam> <reference_fasta> <output_dir> [ --also-vcf <true/false>  --minBaseQ <25> --estimateDispersion --overdispersion <100> --estimateDirichlet --alternative <two.sided> ]"
		exit -1;
	fi

	bed_file=$1; shift # assuming the bed file does not have header lines and just starts with the regions directly, otherwise deepSNV throws and error
	tumor_bam=$1; shift
	normal_bam=$1; shift
    reference_fasta=$1; shift
	out_dir=$1; shift
		
	also_vcf=true
	minBaseQ=25 # default for deepSNV
	estimateDispersion=false
	overdispersion=100 # default for deepSNV model betabinomial
	estimateDirichlet=false
	alternative="two.sided" 
	
	opts=""
	parameter_setting="alternative_${alternative}" # get one string that summarizes all the parameters set - and only specifically mentions them if they are non-default values, except for alternative

	while [ 0 -lt $# ]; do
		if [ $1 == "--also-vcf" ]; then
			shift;
			also_vcf=$1
		elif [ $1 == "--minBaseQ" ]; then
			shift;
			minBaseQ=$1
			if [ $minBaseQ != 25 ]; then
				parameter_setting="${parameter_setting}_minBaseQ${minBaseQ}"
				opts="$opts --minBaseQ $minBaseQ "
			fi
		elif [ $1 == "--estimateDispersion" ]; then
			estimateDispersion=true
			opts="$opts --estimateDispersion "
			parameter_setting="${parameter_setting}_estDispersionTRUE"
		elif [ $1 == "--overdispersion" ]; then
			shift;
			overdispersion=$1
			if [ $overdispersion != 100 ]; then
				parameter_setting="${parameter_setting}_overdisp${overdispersion}"
				opts="$opts --overdispersion $overdispersion "
			fi
		elif [ $1 == "--estimateDirichlet" ]; then
			estimateDirichlet=true
			opts="$opts --estimateDirichlet "
			parameter_setting="${parameter_setting}_estDirichletTRUE"
		elif [ $1 == "--alternative" ]; then
			shift;
			alternative=$1;
			opts="$opts --alternative $alternative "
		fi;
		shift;
	done

	echo "run deepSNV..."
	echo "opts = $opts "

        # check whether the bed file has overlapping regions -  it should not have, because otherwise deepSNV calls variants in these regions several times
        NumOverlaps=`awk 'BEGIN{chr="";end=0}{if(chr==$1 && $2<end){print};chr=$1;end=$3}' ${bed_file} | wc -l `
        if [ $NumOverlaps -gt 0 ]; then
		>&2 echo "ERROR: The bed file $bed_file contains overlapping regions."
		exit -1
        fi

	mkdir -p $out_dir
	
	deepSNV_script="$(dirname ${BASH_SOURCE[0]})/deepSNV_automated.R"
	 
	command="$deepSNV_script $bed_file $tumor_bam $normal_bam $out_dir $opts"
	echo $command
	eval $command

	if $also_vcf; then
		output=$out_dir/`basename ${tumor_bam%.bam}`_vs_`basename ${normal_bam%.bam}`.txt
		fn_vcf=${output%.txt}.vcf
		if [ ! -f $fn_vcf ] ; then 
			command2VCF="/cluster/work/bewi/modules/python/2.7/bin/python $(dirname ${BASH_SOURCE[0]})/deepSNV2VCF.py $output $reference_fasta $fn_vcf --no-LOH --no-indels "
            echo "$command2VCF"
            eval "$command2VCF"
		else
			echo deepSNV vcf file exists: $fn_vcf
		fi
	fi
#}
