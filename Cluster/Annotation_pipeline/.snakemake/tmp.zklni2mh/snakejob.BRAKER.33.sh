#!/bin/sh
# properties = {"params": {"l_mem_free": "4G"}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR0032/"], "rule": "BRAKER", "local": false, "resources": {}, "wildcards": ["BR0032"], "jobid": 33, "threads": 2, "input": ["/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR0032.fasta", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_BR0032.gff"], "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/BR0032/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR0032.fasta /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_BR0032.gff --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/33.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/33.jobfailed"; exit 1)

