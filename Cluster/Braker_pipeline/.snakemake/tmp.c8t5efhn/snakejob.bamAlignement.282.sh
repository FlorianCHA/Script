#!/bin/sh
# properties = {"jobid": 282, "local": false, "log": ["log/bamAlignement_{wildcards.smp}.out"], "rule": "bamAlignement", "params": {}, "input": ["/homedir/charriat/work/Annotation/0_rawdata/rnaseq/", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BF0072.fasta", "/homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt"], "threads": 2, "wildcards": ["BF0072"], "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BF0072/", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BF0072/finalResults/bamList"], "resources": {}, "cluster": {}}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BF0072/ /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/BF0072/finalResults/bamList --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.c8t5efhn /homedir/charriat/work/Annotation/0_rawdata/rnaseq/ /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BF0072.fasta /homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.c8t5efhn/282.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.c8t5efhn/282.jobfailed"; exit 1)

