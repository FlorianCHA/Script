#!/bin/sh
# properties = {"jobid": 4, "threads": 1, "resources": {}, "output": ["/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_AG0004.gff", "/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_AG0004.gff"], "rule": "mergeHint", "cluster": {}, "log": [], "wildcards": ["AG0004"], "params": {}, "input": ["/homedir/charriat/work/Annotation/test/1_hints/RNAseqHints/hints_AG0004.filtered.gff", "/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3"], "local": false}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_AG0004.gff /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_AG0004.gff --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q /homedir/charriat/work/Annotation/test/1_hints/RNAseqHints/hints_AG0004.filtered.gff /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeHint  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q/4.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q/4.jobfailed"; exit 1)

