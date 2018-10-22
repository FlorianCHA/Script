#!/bin/sh
# properties = {"rule": "mergeBam", "local": false, "resources": {}, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1866_scaffold/finalResults/merged_/CH1866_scaffold.accepted_hits.bam"], "cluster": {}, "params": {"l_mem_free": "4G"}, "jobid": 2025, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1866_scaffold/"], "wildcards": ["/CH1866_scaffold"], "log": [], "threads": 1}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1866_scaffold/finalResults/merged_/CH1866_scaffold.accepted_hits.bam --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg /work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1866_scaffold/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeBam  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg/2025.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg/2025.jobfailed"; exit 1)

