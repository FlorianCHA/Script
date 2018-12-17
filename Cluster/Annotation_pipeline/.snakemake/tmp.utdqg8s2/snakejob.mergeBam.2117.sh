#!/bin/sh
# properties = {"resources": {}, "params": {"l_mem_free": "4G"}, "local": false, "rule": "mergeBam", "wildcards": ["CH1164"], "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/0_bamAlignement/CH1164/finalResults/merged_CH1164.accepted_hits.bam"], "cluster": {}, "log": [], "threads": 1, "jobid": 2117, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/0_bamAlignement/CH1164/"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/0_bamAlignement/CH1164/finalResults/merged_CH1164.accepted_hits.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2 /work/gladieux/magMax_project/2_New_Annotation/1_tmp/0_bamAlignement/CH1164/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeBam  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2/2117.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2/2117.jobfailed"; exit 1)
