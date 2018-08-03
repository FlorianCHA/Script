#!/bin/sh
# properties = {"local": false, "rule": "exonerate", "jobid": 194, "threads": 2, "cluster": {}, "input": ["SupplementaryFile/70-15_proteinAndAVR.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/TN0057.fasta"], "params": {}, "wildcards": ["TN0057"], "log": [], "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_TN0057.hints.gff3", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_TN0057.gff3"], "resources": {}}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_TN0057.hints.gff3 /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_TN0057.gff3 --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.7nvqeb7_ SupplementaryFile/70-15_proteinAndAVR.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/TN0057.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.7nvqeb7_/194.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.7nvqeb7_/194.jobfailed"; exit 1)

