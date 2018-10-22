#!/bin/sh
# properties = {"input": ["/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR29.fasta", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_BR29.gff"], "log": [], "local": false, "threads": 2, "jobid": 32, "rule": "BRAKER", "resources": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR29/"], "params": {}, "wildcards": ["BR29"], "cluster": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR29/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR29.fasta /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_BR29.gff --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl/32.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl/32.jobfailed"; exit 1)

