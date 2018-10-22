#!/bin/sh
# properties = {"cluster": {}, "rule": "sortBam", "log": [], "wildcards": ["BF0072"], "params": {"l_mem_free": "40G"}, "threads": 1, "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits_sort.bam"], "jobid": 24, "local": false, "resources": {}, "input": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits.bam"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits_sort.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits.bam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules sortBam  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy/24.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.gb1f3gxy/24.jobfailed"; exit 1)

