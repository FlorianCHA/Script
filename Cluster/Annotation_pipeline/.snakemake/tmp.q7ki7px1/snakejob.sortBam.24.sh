#!/bin/sh
# properties = {"jobid": 24, "input": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits.bam"], "threads": 1, "wildcards": ["BN0019"], "local": false, "params": {"l_mem_free": "40G"}, "resources": {}, "log": [], "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits_sort.bam"], "rule": "sortBam", "cluster": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits_sort.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1 /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits.bam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules sortBam  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1/24.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1/24.jobfailed"; exit 1)

