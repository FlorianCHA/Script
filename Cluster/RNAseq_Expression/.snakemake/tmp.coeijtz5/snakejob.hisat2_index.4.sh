#!/bin/sh
# properties = {"rule": "hisat2_index", "local": false, "resources": {}, "input": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/Genome/Y34_merge.gff3"], "jobid": 4, "wildcards": ["S117"], "cluster": {}, "params": {"l_mem_free": "4G"}, "threads": 1, "output": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/Ref/S117"], "log": []}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/Ref/S117 --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.coeijtz5 /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/Genome/Y34_merge.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules hisat2_index  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.coeijtz5/4.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.coeijtz5/4.jobfailed"; exit 1)

