#!/bin/sh
# properties = {"cluster": {}, "rule": "bam2hints", "log": [], "wildcards": ["BN0019"], "params": {}, "threads": 2, "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/hints_BN0019.raw.bam"], "jobid": 20, "local": false, "resources": {}, "input": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits_sort.bam"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/hints_BN0019.raw.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BN0019/finalResults/merged_BN0019.accepted_hits_sort.bam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bam2hints  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy/20.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy/20.jobfailed"; exit 1)

