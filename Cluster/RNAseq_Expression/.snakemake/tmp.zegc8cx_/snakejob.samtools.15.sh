#!/bin/sh
# properties = {"cluster": {}, "params": {"l_mem_free": "4G"}, "log": [], "output": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001.bam", "/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001_sort.bam"], "jobid": 15, "threads": 1, "local": false, "resources": {}, "wildcards": ["CL367_c_S16_L005_R1_001"], "rule": "samtools", "input": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/sam_file/CL367_c_S16_L005_R1_001.sam"]}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001.bam /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001_sort.bam --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_ /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/sam_file/CL367_c_S16_L005_R1_001.sam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules samtools  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_/15.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_/15.jobfailed"; exit 1)

