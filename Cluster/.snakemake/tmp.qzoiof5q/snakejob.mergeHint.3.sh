#!/bin/sh
# properties = {"jobid": 3, "threads": 1, "resources": {}, "output": ["/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_CH1857.gff", "/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_CH1857.gff"], "rule": "mergeHint", "cluster": {}, "log": [], "wildcards": ["CH1857"], "params": {}, "input": ["/homedir/charriat/work/Annotation/test/1_hints/RNAseqHints/hints_CH1857.filtered.gff", "/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.hints.gff3"], "local": false}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_CH1857.gff /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_CH1857.gff --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q /homedir/charriat/work/Annotation/test/1_hints/RNAseqHints/hints_CH1857.filtered.gff /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.hints.gff3 --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules mergeHint  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q/3.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.qzoiof5q/3.jobfailed"; exit 1)

