#!/bin/sh
# properties = {"resources": {}, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints.intron_B71.gff", "/homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/B71.fasta"], "cluster": {}, "threads": 2, "jobid": 895, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/2_Braker/B71/"], "rule": "BRAKER", "local": false, "params": {"l_mem_free": "10G"}, "wildcards": ["B71"], "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/2_Braker/B71/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints.intron_B71.gff /homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/B71.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr/895.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr/895.jobfailed"; exit 1)

