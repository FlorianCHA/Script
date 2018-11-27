#!/bin/sh
# properties = {"resources": {}, "input": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints_PH0014-rn.gff", "/homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/PH0014-rn.fasta"], "cluster": {}, "threads": 2, "jobid": 738, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/3_augustus/PH0014-rn.gff3"], "rule": "augustus", "local": false, "params": {"l_mem_free": "4G"}, "wildcards": ["PH0014-rn"], "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/3_augustus/PH0014-rn.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/MergeHints/RNAseq_protein.hints_PH0014-rn.gff /homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/PH0014-rn.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules augustus  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr/738.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp.pxe06fkr/738.jobfailed"; exit 1)

