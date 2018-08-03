#!/bin/sh
# properties = {"resources": {}, "threads": 2, "jobid": 327, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BN0202/"], "rule": "bamAlignement", "params": {}, "wildcards": ["BN0202"], "log": ["log/bamAlignement_{wildcards.smp}.out"], "cluster": {}, "input": ["SupplementaryFile/tophatMapping.config.txt", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BN0202.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq"], "local": false}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BN0202/ --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.8bprgh_1 SupplementaryFile/tophatMapping.config.txt /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BN0202.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.8bprgh_1/327.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.8bprgh_1/327.jobfailed"; exit 1)

