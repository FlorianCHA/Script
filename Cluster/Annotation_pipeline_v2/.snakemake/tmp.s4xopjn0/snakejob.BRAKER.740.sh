#!/bin/sh
# properties = {"threads": 2, "local": false, "input": ["/work/gladieux/magMax_project/2_New_Annotation/0_rawdata/IN0116.fasta", "/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints.intron_IN0116.gff"], "log": [], "cluster": {}, "rule": "BRAKER", "wildcards": ["IN0116"], "jobid": 740, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/2_Braker/IN0116/"], "params": {"l_mem_free": "10G", "species": "magnaporthe_oryzae"}, "resources": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline_v2 && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/2_Braker/IN0116/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline_v2/BRAKER_pipeline2.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0 /work/gladieux/magMax_project/2_New_Annotation/0_rawdata/IN0116.fasta /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints.intron_IN0116.gff --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/740.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/740.jobfailed"; exit 1)

