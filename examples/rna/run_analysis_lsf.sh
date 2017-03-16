NAKEMAKE=$1

module load java; module load R; module load python/2.7.6; bsub -n 1 -W 235 -R "rusage[mem=2000]" -e rna.err -o rna.out $SNAKEMAKE -s rna.snake --configfile rna_config.json --cluster "bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 3000 -p
