#!/bin/sh
# properties = {"jobid": 18, "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/FR13/", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/FR13.gff3"], "wildcards": ["FR13"], "params": {"l_mem_free": "4G"}, "resources": {}, "log": [], "cluster": {}, "local": false, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/FR13_merge.gff3"], "rule": "mergeGFF", "threads": 1}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/4_mergeGFF/FR13_merge.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/FR13/ /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/FR13.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeGFF  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/18.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/18.jobfailed"; exit 1)

