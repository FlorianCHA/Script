#!/bin/sh
# properties = {"output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/VT0027/"], "wildcards": ["VT0027"], "jobid": 299, "local": false, "resources": {}, "rule": "bamAlignement", "params": {}, "cluster": {}, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq", "SupplementaryFile/tophatMapping.config.txt", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/VT0027.fasta"], "log": ["log/bamAlignement_{wildcards.smp}.out"], "threads": 2}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/VT0027/ --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.ibc21s3b /work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq SupplementaryFile/tophatMapping.config.txt /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/VT0027.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.ibc21s3b/299.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.ibc21s3b/299.jobfailed"; exit 1)

