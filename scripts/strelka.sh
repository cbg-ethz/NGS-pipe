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
configure_strelka_workflow=${10}
gatk=${11}  # gatk + gatk_key
#sourcePerl=${12}

#echo $1,$2,$3,$4,$5,$6,$7,$8,$9,${10},${11},${12}

mkdir -p $project_dir/temp_$label
tmp_dir=$project_dir/temp_$label
out_dir=$project_dir

##### Run strelka

$configure_strelka_workflow \
--tumor $tumor_bam \
--normal $normal_bam \
--ref $reference \
--config $strelka_config \
--output-dir $tmp_dir/strelka

make -C $tmp_dir/strelka

### SNVs

#echo "GATK CALL:" $gatk
$gatk \
-T SelectVariants \
-R $reference \
--variant $tmp_dir/strelka/results/all.somatic.snvs.vcf \
-L $intervals \
--interval_padding $interval_padding \
--out $out_dir/$label.snv.vcf

### Indels

$gatk \
-T SelectVariants \
-R $reference \
--variant $tmp_dir/strelka/results/all.somatic.indels.vcf \
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
