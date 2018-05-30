#!/bin/sh
# properties = {"jobid": 14, "local": false, "cluster": {}, "params": {}, "wildcards": ["AG0004"], "output": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004/merged_AG0004.accepted_hits.bam"], "rule": "mergeBam", "input": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004/"], "threads": 2, "resources": {}, "log": []}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004/merged_AG0004.accepted_hits.bam --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.m3g2lcuh /homedir/charriat/work/Annotation/test/0_bamAlignement/AG0004/ --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeBam  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.m3g2lcuh/14.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.m3g2lcuh/14.jobfailed"; exit 1)

