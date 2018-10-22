#!/bin/sh
# properties = {"cluster": {}, "wildcards": ["TH0016"], "resources": {}, "rule": "augustus", "local": false, "params": {"l_mem_free": "4G"}, "jobid": 33, "input": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints_TH0016.gff", "/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta"], "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/TH0016.gff3"], "threads": 2, "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/3_augustus/TH0016.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ozzl3qpk /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints_TH0016.gff /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules augustus  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ozzl3qpk/33.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.ozzl3qpk/33.jobfailed"; exit 1)

