#!/bin/sh
# properties = {"params": {"l_mem_free": "4G"}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/BR0032_merge.gff3"], "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/BR0032.gff3", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR0032/"], "rule": "mergeGFF", "threads": 1, "wildcards": ["BR0032"], "jobid": 17, "resources": {}, "local": false, "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/BR0032_merge.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.7xjrwakn /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/BR0032.gff3 /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR0032/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeGFF  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.7xjrwakn/17.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.7xjrwakn/17.jobfailed"; exit 1)

