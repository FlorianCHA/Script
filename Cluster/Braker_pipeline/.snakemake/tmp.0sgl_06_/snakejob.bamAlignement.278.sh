#!/bin/sh
# properties = {"resources": {}, "local": false, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/IN0017.fasta", "/homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt", "/homedir/charriat/work/Annotation/0_rawdata/rnaseq/"], "wildcards": ["IN0017"], "params": {}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/IN0017"], "threads": 2, "rule": "bamAlignement", "log": ["log/bamAlignement_{wildcards.smp}.out"], "jobid": 278}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/0_bamAlignement/IN0017 --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_ /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/IN0017.fasta /homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt /homedir/charriat/work/Annotation/0_rawdata/rnaseq/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_/278.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_/278.jobfailed"; exit 1)

