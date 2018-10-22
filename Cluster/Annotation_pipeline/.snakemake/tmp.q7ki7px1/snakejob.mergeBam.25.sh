#!/bin/sh
# properties = {"jobid": 25, "input": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/"], "threads": 1, "wildcards": ["BF0072"], "local": false, "params": {}, "resources": {}, "log": [], "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits.bam"], "rule": "mergeBam", "cluster": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/finalResults/merged_BF0072.accepted_hits.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1 /work/gladieux/magMax_project/2_Annotation/6_test/0_bamAlignement/BF0072/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeBam  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1/25.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.q7ki7px1/25.jobfailed"; exit 1)

