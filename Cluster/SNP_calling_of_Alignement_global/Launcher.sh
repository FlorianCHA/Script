#!/bin/bash
module purge
module load bioinfo/snakemake/3.13.3
snakemake -s BRAKER_pipeline.snake --jobs 100 --cluster "qsub -q long.q -cwd -V -pe parallel_smp {threads} -l mem_free={params.l_mem_free}"
