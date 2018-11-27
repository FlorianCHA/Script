#!/bin/sh
# properties = {"jobid": 1150, "log": [], "input": ["SupplementaryFile/OG.fasta", "/homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/GFSI1-7-2.fasta"], "rule": "exonerate", "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonarate_GFSI1-7-2.gff3", "/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_GFSI1-7-2.hints.gff3"], "cluster": {}, "resources": {}, "wildcards": ["GFSI1-7-2"], "local": false, "threads": 2, "params": {"l_mem_free": "4G"}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonarate_GFSI1-7-2.gff3 /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_GFSI1-7-2.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.1vrdbnwr SupplementaryFile/OG.fasta /homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/GFSI1-7-2.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.1vrdbnwr/1150.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.1vrdbnwr/1150.jobfailed"; exit 1)

