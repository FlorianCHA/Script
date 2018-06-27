#!/bin/sh
# properties = {"output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BR29/"], "wildcards": ["BR29"], "rule": "bamAlignement", "resources": {}, "threads": 2, "jobid": 269, "cluster": {}, "log": ["log/bamAlignement_{wildcards.smp}.out"], "local": false, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BR29.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq", "SupplementaryFile/tophatMapping.config.txt"], "params": {}}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BR29/ --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BR29.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq SupplementaryFile/tophatMapping.config.txt --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc/269.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc/269.jobfailed"; exit 1)

