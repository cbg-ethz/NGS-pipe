#!/bin/bash

SNAKEMAKE=$1

module load java; bsub -J SRA053195 -n 1 -W 119:59 -R "rusage[mem=20000]" -e SRP051153.err -o SRP051153.out $SNAKEMAKE -s SRP051153.snake --configfile config.json --cluster "bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 3000 -p
