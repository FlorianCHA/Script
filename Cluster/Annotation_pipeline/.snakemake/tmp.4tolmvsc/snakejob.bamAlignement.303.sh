#!/bin/sh
# properties = {"output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/"], "wildcards": ["CH1205-canu2"], "rule": "bamAlignement", "resources": {}, "threads": 2, "jobid": 303, "cluster": {}, "log": ["log/bamAlignement_{wildcards.smp}.out"], "local": false, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/CH1205-canu2.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq", "SupplementaryFile/tophatMapping.config.txt"], "params": {}}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/ --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/CH1205-canu2.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq SupplementaryFile/tophatMapping.config.txt --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc/303.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.4tolmvsc/303.jobfailed"; exit 1)

