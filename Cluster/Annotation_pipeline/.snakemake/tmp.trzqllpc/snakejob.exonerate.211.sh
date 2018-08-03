#!/bin/sh
# properties = {"resources": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_IN0059.gff3", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_IN0059.hints.gff3"], "threads": 4, "jobid": 211, "local": false, "log": [], "cluster": {}, "wildcards": ["IN0059"], "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/IN0059.fasta", "SupplementaryFile/70-15_annotated_protein.fa"], "rule": "exonerate", "params": {}}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_IN0059.gff3 /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_IN0059.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.trzqllpc /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/IN0059.fasta SupplementaryFile/70-15_annotated_protein.fa --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.trzqllpc/211.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.trzqllpc/211.jobfailed"; exit 1)

