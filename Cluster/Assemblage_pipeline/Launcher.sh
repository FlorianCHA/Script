#!/bin/bash
module unload system/python/2.7.9
module load bioinfo/snakemake/3.13.3
snakemake -s ABySS_pipeline.snake --jobs 30 --cluster "qsub -q normal.q -cwd -V -pe parallel_smp {threads} -l mem_free=40G"
