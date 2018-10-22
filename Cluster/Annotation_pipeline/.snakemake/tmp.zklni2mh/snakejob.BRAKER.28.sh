#!/bin/sh
# properties = {"params": {"l_mem_free": "4G"}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/TH0016/"], "rule": "BRAKER", "local": false, "resources": {}, "wildcards": ["TH0016"], "jobid": 28, "threads": 2, "input": ["/home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta", "/work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_TH0016.gff"], "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/2_Braker/TH0016/ --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh /home/gladieux/work/magMax_project/2_Annotation/0_rawdata/assembly_mBio/GEMO/TH0016.fasta /work/gladieux/magMax_project/2_Annotation/6_GEMO_TEST/1_hints/MergeHints/RNAseq_protein.hints.intron_TH0016.gff --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/28.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.zklni2mh/28.jobfailed"; exit 1)

