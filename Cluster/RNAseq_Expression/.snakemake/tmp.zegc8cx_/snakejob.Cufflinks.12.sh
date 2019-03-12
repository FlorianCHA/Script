#!/bin/sh
# properties = {"cluster": {}, "params": {"l_mem_free": "4G"}, "log": [], "output": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/cufflinks/CL367_d_S17_L005_R1_001"], "jobid": 12, "threads": 1, "local": false, "resources": {}, "wildcards": ["CL367_d_S17_L005_R1_001"], "rule": "Cufflinks", "input": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_d_S17_L005_R1_001_sort.bam", "/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/genome/CH1908_merge.gff3"]}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/cufflinks/CL367_d_S17_L005_R1_001 --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_ /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_d_S17_L005_R1_001_sort.bam /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/genome/CH1908_merge.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules Cufflinks  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_/12.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.zegc8cx_/12.jobfailed"; exit 1)

