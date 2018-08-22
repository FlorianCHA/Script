#!/bin/sh
# properties = {"wildcards": ["AG0004"], "threads": 2, "jobid": 2, "resources": {}, "local": false, "params": {}, "input": ["/homedir/charriat/work/Annotation/0_rawdata/rnaseq/", "/homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt", "/homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta"], "rule": "bamAlignement", "output": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004"], "log": [], "cluster": {}}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004 --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_ /homedir/charriat/work/Annotation/0_rawdata/rnaseq/ /homedir/charriat/work/Annotation/1_tmp/ToogleConfig/tophatMapping.config.txt /homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_/2.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.yp3sllu_/2.jobfailed"; exit 1)

