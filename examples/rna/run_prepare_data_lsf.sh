SNAKEMAKE=$1

bsub -J prepare -n 1 -W 235 -R "rusage[mem=2000]" -e prepare_data.err -o prepare_data.out $SNAKEMAKE -s prepare_data.snake --cluster "bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 100 -p -k

