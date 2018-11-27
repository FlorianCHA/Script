#!/bin/sh
# properties = {"resources": {}, "params": {"l_mem_free": "4G"}, "local": false, "rule": "renameFasta", "wildcards": ["P29"], "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/5_final_data/P29/P29_protein.fasta"], "cluster": {}, "log": [], "threads": 1, "jobid": 98, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/5_final_data/P29/"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/5_final_data/P29/P29_protein.fasta --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2 /work/gladieux/magMax_project/2_New_Annotation/1_tmp/5_final_data/P29/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules renameFasta  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2/98.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.utdqg8s2/98.jobfailed"; exit 1)

