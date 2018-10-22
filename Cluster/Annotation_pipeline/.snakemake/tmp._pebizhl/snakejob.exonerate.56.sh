#!/bin/sh
# properties = {"input": ["SupplementaryFile/70-15_proteinAndAVR.fasta", "/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/US0071.fasta"], "log": [], "local": false, "threads": 2, "jobid": 56, "rule": "exonerate", "resources": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_US0071.gff3", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_US0071.hints.gff3"], "params": {"l_mem_free": "4G"}, "wildcards": ["US0071"], "cluster": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_US0071.gff3 /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/ProtHints/exonarate_US0071.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl SupplementaryFile/70-15_proteinAndAVR.fasta /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/US0071.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl/56.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._pebizhl/56.jobfailed"; exit 1)

