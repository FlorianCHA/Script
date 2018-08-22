#!/bin/sh
# properties = {"local": false, "cluster": {}, "rule": "tblastn", "wildcards": ["OG0000039"], "threads": 1, "params": {"name": "All_genome_isolat"}, "log": [], "jobid": 5321, "input": ["/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/1_test_snakemake/tmp/0_OG-fasta/OG0000039.fasta", "/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/0_DB/"], "resources": {}, "output": ["/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/1_test_snakemake/tmp/1_tblastn-result/OG0000039_blast.txt"]}
cd /work/gladieux/Script/Cluster/Correction_Orthofinder && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/1_test_snakemake/tmp/1_tblastn-result/OG0000039_blast.txt --snakefile /work/gladieux/Script/Cluster/Correction_Orthofinder/correctionOrthofinder.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Correction_Orthofinder/.snakemake/tmp.c5bfho4s /homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/1_test_snakemake/tmp/0_OG-fasta/OG0000039.fasta /homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/0_DB/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules tblastn  && touch "/work/gladieux/Script/Cluster/Correction_Orthofinder/.snakemake/tmp.c5bfho4s/5321.jobfinished" || (touch "/work/gladieux/Script/Cluster/Correction_Orthofinder/.snakemake/tmp.c5bfho4s/5321.jobfailed"; exit 1)

