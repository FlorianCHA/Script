#!/bin/bash
module unload system/python/2.7.9
module load bioinfo/snakemake/3.13.3
snakemake -s BRAKER_pipeline.snake --jobs 25 --cluster "qsub -q long.q -l mem_free=30G -cwd -V -pe parallel_smp {threads}"
