#!/bin/sh
# properties = {"resources": {}, "log": [], "cluster": {}, "jobid": 12, "wildcards": ["CH1857"], "input": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857/merged_CH1857.accepted_hits.bam"], "threads": 2, "output": ["/homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857/merged_CH1857.accepted_hits_sort.bam"], "rule": "sortBam", "params": {}, "local": false}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857/merged_CH1857.accepted_hits_sort.bam --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.egdbz33l /homedir/charriat/work/Annotation/test/0_bamAlignement/CH1857/merged_CH1857.accepted_hits.bam --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules sortBam  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.egdbz33l/12.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.egdbz33l/12.jobfailed"; exit 1)

