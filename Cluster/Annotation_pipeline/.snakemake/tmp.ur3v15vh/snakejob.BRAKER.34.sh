#!/bin/sh
# properties = {"jobid": 34, "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_US0071.gff", "/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/US0071.fasta"], "wildcards": ["US0071"], "params": {"l_mem_free": "4G"}, "resources": {}, "log": [], "cluster": {}, "local": false, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/US0071/"], "rule": "BRAKER", "threads": 2}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/US0071/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_US0071.gff /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/US0071.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/34.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/34.jobfailed"; exit 1)

