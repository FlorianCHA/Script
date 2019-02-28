#!/bin/sh
# properties = {"local": false, "output": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/cufflinks/S120"], "rule": "Cufflinks", "wildcards": ["S120"], "log": [], "resources": {}, "jobid": 6, "input": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S120_sort.bam", "/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/Genome/Y34_merge.gff3"], "params": {"l_mem_free": "4G"}, "threads": 1, "cluster": {}}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/cufflinks/S120 --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.369midka /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S120_sort.bam /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/Genome/Y34_merge.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules Cufflinks  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.369midka/6.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.369midka/6.jobfailed"; exit 1)

