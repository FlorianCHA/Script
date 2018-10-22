#!/bin/sh
# properties = {"params": {}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/BR29_merge.gff3"], "rule": "mergeGFF", "local": false, "resources": {}, "wildcards": ["BR29"], "jobid": 22, "threads": 1, "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR29/", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/BR29.gff3"], "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/BR29_merge.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR29/ /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/BR29.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeGFF  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/22.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/22.jobfailed"; exit 1)

