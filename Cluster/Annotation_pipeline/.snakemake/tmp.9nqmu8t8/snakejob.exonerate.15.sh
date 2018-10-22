#!/bin/sh
# properties = {"log": [], "rule": "exonerate", "wildcards": ["BN0019"], "params": {}, "jobid": 15, "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.hints.gff3", "/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.gff3"], "threads": 2, "cluster": {}, "local": false, "resources": {}, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BN0019.fasta", "SupplementaryFile/70-15_proteinAndAVR.fasta"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.hints.gff3 /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.9nqmu8t8 /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BN0019.fasta SupplementaryFile/70-15_proteinAndAVR.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.9nqmu8t8/15.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.9nqmu8t8/15.jobfailed"; exit 1)

