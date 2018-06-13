#!/bin/sh
# properties = {"rule": "bamAlignement", "cluster": {}, "log": ["log/bamAlignement_{wildcards.smp}.out"], "input": ["/homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt", "/homedir/charriat/work/Annotation/0_rawdata/rnaseq/", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/CH1205-canu2.fasta"], "params": {}, "wildcards": ["CH1205-canu2"], "threads": 2, "resources": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/finalResults/bamList"], "jobid": 274, "local": false}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/ /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/CH1205-canu2/finalResults/bamList --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07 /homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt /homedir/charriat/work/Annotation/0_rawdata/rnaseq/ /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/CH1205-canu2.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07/274.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07/274.jobfailed"; exit 1)

