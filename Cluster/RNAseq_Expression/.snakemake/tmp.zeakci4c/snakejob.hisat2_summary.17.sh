#!/bin/sh
# properties = {"resources": {}, "threads": 1, "cluster": {}, "output": ["/work/gladieux/magMax_project/5_RNAseq/BR32/result/sam_file/46_ATGTCA_L005.sam", "/work/gladieux/magMax_project/5_RNAseq/BR32/result/sam_file/46_ATGTCA_L005_cufflinks.txt"], "log": [], "rule": "hisat2_summary", "local": false, "wildcards": ["46_ATGTCA_L005"], "params": {"l_mem_free": "4G"}, "jobid": 17, "input": ["/work/gladieux/magMax_project/5_RNAseq/BR32/RNAseq/46_ATGTCA_L005_R1.fastq.gz", "/work/gladieux/magMax_project/5_RNAseq/BR32/Genome/BR0032", "/work/gladieux/magMax_project/5_RNAseq/BR32/RNAseq/46_ATGTCA_L005_R2.fastq.gz"]}
cd /gs7k1/home/gladieux/data_flo/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/5_RNAseq/BR32/result/sam_file/46_ATGTCA_L005.sam /work/gladieux/magMax_project/5_RNAseq/BR32/result/sam_file/46_ATGTCA_L005_cufflinks.txt --snakefile /gs7k1/home/gladieux/data_flo/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/gladieux/data_flo/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zeakci4c /work/gladieux/magMax_project/5_RNAseq/BR32/RNAseq/46_ATGTCA_L005_R1.fastq.gz /work/gladieux/magMax_project/5_RNAseq/BR32/Genome/BR0032 /work/gladieux/magMax_project/5_RNAseq/BR32/RNAseq/46_ATGTCA_L005_R2.fastq.gz --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules hisat2_summary  && touch "/gs7k1/home/gladieux/data_flo/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zeakci4c/17.jobfinished" || (touch "/gs7k1/home/gladieux/data_flo/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zeakci4c/17.jobfailed"; exit 1)

