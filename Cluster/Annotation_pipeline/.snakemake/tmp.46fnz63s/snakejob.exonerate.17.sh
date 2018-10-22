#!/bin/sh
# properties = {"log": [], "input": ["SupplementaryFile/70-15_proteinAndAVR.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BF0072.fasta"], "jobid": 17, "params": {}, "local": false, "threads": 2, "resources": {}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BF0072.hints.gff3", "/work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BF0072.gff3"], "rule": "exonerate", "wildcards": ["BF0072"]}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BF0072.hints.gff3 /work/gladieux/magMax_project/2_Annotation/6_test/1_hints/ProtHints/exonarate_BF0072.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.46fnz63s SupplementaryFile/70-15_proteinAndAVR.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/AssemblyToulouse/test/BF0072.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.46fnz63s/17.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.46fnz63s/17.jobfailed"; exit 1)

