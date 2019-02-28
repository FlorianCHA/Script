#!/bin/sh
# properties = {"log": [], "wildcards": ["S117"], "jobid": 9, "resources": {}, "rule": "samtools", "cluster": {}, "input": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/sam_file/S117.sam"], "threads": 1, "params": {"l_mem_free": "4G"}, "local": false, "output": ["/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S117.bam", "/homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S117_sort.bam"]}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S117.bam /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/bam_file/S117_sort.bam --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.9imu1y4_ /homedir/gladieux/work/magMax_project/5_RNAseq/1_Alignement/test/0_rawdata/result/sam_file/S117.sam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules samtools  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.9imu1y4_/9.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.9imu1y4_/9.jobfailed"; exit 1)

