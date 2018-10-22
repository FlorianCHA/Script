#!/bin/sh
# properties = {"jobid": 36, "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_TH0016.gff", "/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta"], "wildcards": ["TH0016"], "params": {"l_mem_free": "4G"}, "resources": {}, "log": [], "cluster": {}, "local": false, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/TH0016/"], "rule": "BRAKER", "threads": 2}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/TH0016/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_TH0016.gff /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/36.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ur3v15vh/36.jobfailed"; exit 1)

