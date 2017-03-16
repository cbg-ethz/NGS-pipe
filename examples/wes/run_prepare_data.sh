#!/bin/bash

SNAKEMAKE=$1

bsub -J prepare -n 1 -W 23:59 -R "rusage[mem=20000]" -e prepare_data.err -o prepare_data.out $SNAKEMAKE -s prepare_data.snake --cluster "bsub -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 100 -p -k
