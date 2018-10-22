#!/bin/sh
# properties = {"log": [], "input": ["/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR29.fasta", "SupplementaryFile/70-15_proteinAndAVR.fasta"], "rule": "exonerate", "resources": {}, "wildcards": ["BR29"], "local": false, "jobid": 50, "params": {"l_mem_free": "4G"}, "threads": 2, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_BR29.gff3", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_BR29.hints.gff3"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_BR29.gff3 /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_BR29.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.56l5520t /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/BR29.fasta SupplementaryFile/70-15_proteinAndAVR.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.56l5520t/50.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.56l5520t/50.jobfailed"; exit 1)

