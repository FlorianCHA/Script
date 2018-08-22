#!/bin/sh
# properties = {"wildcards": ["CH1857"], "threads": 2, "jobid": 1, "resources": {}, "local": false, "params": {}, "input": ["/homedir/charriat/work/Annotation/0_rawdata/rnaseq/", "/homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt", "/homedir/charriat/work/Annotation/test_input/CH1857_scaffold.fasta"], "rule": "bamAlignement", "output": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857"], "log": [], "cluster": {}}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857 --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_ /homedir/charriat/work/Annotation/0_rawdata/rnaseq/ /homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt /homedir/charriat/work/Annotation/test_input/CH1857_scaffold.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_/1.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_/1.jobfailed"; exit 1)

