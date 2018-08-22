#!/bin/bash
module unload system/python/2.7.9
module load bioinfo/snakemake/3.13.3
snakemake -s ABySS_pipeline.snake --dag | dot -Tpdf > dag.pdf
