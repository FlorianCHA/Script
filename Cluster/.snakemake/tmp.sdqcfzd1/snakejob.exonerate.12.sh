#!/bin/sh
# properties = {"output": ["/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.gff3", "/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3"], "jobid": 12, "params": {}, "local": false, "rule": "exonerate", "cluster": {}, "resources": {}, "wildcards": ["AG0004"], "log": [], "threads": 4, "input": ["/homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa", "/homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta"]}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.gff3 /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.sdqcfzd1 /homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa /homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.sdqcfzd1/12.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.sdqcfzd1/12.jobfailed"; exit 1)

