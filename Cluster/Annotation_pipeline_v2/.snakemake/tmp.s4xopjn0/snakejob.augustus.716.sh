#!/bin/sh
# properties = {"threads": 2, "local": false, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints_VT0027.gff", "/work/gladieux/magMax_project/2_New_Annotation/0_rawdata/VT0027.fasta"], "log": [], "cluster": {}, "rule": "augustus", "wildcards": ["VT0027"], "jobid": 716, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/3_augustus/VT0027.gff3"], "params": {"l_mem_free": "4G", "species": "magnaporthe_oryzae"}, "resources": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline_v2 && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/3_augustus/VT0027.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline_v2/BRAKER_pipeline2.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0 /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints_VT0027.gff /work/gladieux/magMax_project/2_New_Annotation/0_rawdata/VT0027.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules augustus  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/716.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/716.jobfailed"; exit 1)

