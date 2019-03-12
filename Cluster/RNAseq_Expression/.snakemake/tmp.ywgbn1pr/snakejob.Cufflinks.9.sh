#!/bin/sh
# properties = {"wildcards": ["CL367_c_S16_L005_R1_001"], "log": [], "cluster": {}, "jobid": 9, "input": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001_sort.bam", "/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/genome/CH1908_merge.gff3"], "params": {"l_mem_free": "4G"}, "threads": 1, "resources": {}, "local": false, "output": ["/work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/cufflinks/CL367_c_S16_L005_R1_001"], "rule": "Cufflinks"}
cd /work/gladieux/Script/Cluster/RNAseq_Expression && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/cufflinks/CL367_c_S16_L005_R1_001 --snakefile /work/gladieux/Script/Cluster/RNAseq_Expression/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.ywgbn1pr /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/result/bam_file/CL367_c_S16_L005_R1_001_sort.bam /work/gladieux/magMax_project/5_RNAseq/1_Alignement/CL367/genome/CH1908_merge.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules Cufflinks  && touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.ywgbn1pr/9.jobfinished" || (touch "/work/gladieux/Script/Cluster/RNAseq_Expression/.snakemake/tmp.ywgbn1pr/9.jobfailed"; exit 1)

