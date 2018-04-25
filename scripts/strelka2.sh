#!/bin/bash

project_dir=$1
reference=$2
dbsnp=$3
intervals=$4
interval_padding=$5

strelka_config=$6

tumor_bam=$7
normal_bam=$8

label=$9
configure_strelka_workflow_1=${10}
gatk=${11}  # gatk + gatk_key
#sourcePerl=${12}

#echo $1,$2,$3,$4,$5,$6,$7,$8,$9,${10},${11},${12}

mkdir -p $project_dir/temp_$label
tmp_dir=$project_dir/temp_$label
out_dir=$project_dir

#source $sourcePerl
#module load vcp
#source /cluster/apps/vcp/perl5/etc/bashrc
#configure_strelka_workflow=/cluster/apps/vcp/strelka/1.0.14/x86_64/bin/configureStrelkaWorkflow.pl

#configure_strelka_workflow=$(readlink $configure_strelka_workflow_1)
configure_strelka_workflow=$configure_strelka_workflow_1
##### Run strelka
echo "$configure_strelka_workflow --tumorBam $tumor_bam --normalBam $normal_bam --ref $reference --config $strelka_config --exome --runDir $tmp_dir/strelka"
$configure_strelka_workflow --tumorBam $tumor_bam --normalBam $normal_bam --ref $reference --config $strelka_config --exome --runDir $tmp_dir/strelka

#make -C $tmp_dir/strelka

python2 $tmp_dir/strelka/runWorkflow.py --mode local

### SNVs

#echo "GATK CALL:" $gatk
$gatk \
-T SelectVariants \
-R $reference \
--variant $tmp_dir/strelka/results/variants/somatic.snvs.vcf.gz \
-L $intervals \
--interval_padding $interval_padding \
--out $out_dir/$label.snv.vcf

### Indels

$gatk \
-T SelectVariants \
-R $reference \
--variant $tmp_dir/strelka/results/variants/somatic.indels.vcf.gz \
-L $intervals \
--interval_padding $interval_padding \
--out $out_dir/$label.indel.vcf

##### Combine indel and snv

$gatk \
-T CombineVariants \
-R $reference \
--variant:strelkaSNV $out_dir/$label.snv.vcf \
--variant:strelkaIndel $out_dir/$label.indel.vcf \
-o $out_dir/$label.vcf \
-genotypeMergeOptions PRIORITIZE \
-priority strelkaSNV,strelkaIndel
