#!/bin/sh
# properties = {"cluster": {}, "params": {"l_mem_free": "4G"}, "resources": {}, "threads": 1, "output": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/IN0059/Genome/IN0059"], "log": [], "wildcards": [], "input": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/IN0059/Genome/IN0059.fasta"], "rule": "hisat2_index", "local": false, "jobid": 14}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/5_RNAseq/1_Alignement/IN0059/Genome/IN0059 --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.nkf2r3rk /work/gladieux/magMax_project/5_RNAseq/1_Alignement/IN0059/Genome/IN0059.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules hisat2_index  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.nkf2r3rk/14.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.nkf2r3rk/14.jobfailed"; exit 1)

