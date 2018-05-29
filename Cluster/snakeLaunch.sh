module unload system/python/2.7.9
module load bioinfo/snakemake/3.13.3
snakemake -s BRAKER_pipeline.snake --jobs 50  --cluster "qsub -q long.q -cwd -V -l mem_free=30G -pe parallel_smp {threads}"
