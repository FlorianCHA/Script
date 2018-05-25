snakemake -s BRAKER_pipeline.snake --jobs 50 --drmaa-log-dir ./trash --cluster "qsub -q long.q -cwd -V -l mem_free=30G -pe parallel_smp {threads}"
