#!/bin/sh
# properties = {"rule": "exonerate", "params": {}, "wildcards": ["BN0019"], "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.gff3", "/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.hints.gff3"], "cluster": {}, "threads": 2, "input": ["SupplementaryFile/70-15_proteinAndAVR.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BN0019.fasta"], "log": [], "resources": {}, "local": false, "jobid": 18}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.gff3 /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BN0019.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.qpmbs7gd SupplementaryFile/70-15_proteinAndAVR.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BN0019.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.qpmbs7gd/18.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.qpmbs7gd/18.jobfailed"; exit 1)

