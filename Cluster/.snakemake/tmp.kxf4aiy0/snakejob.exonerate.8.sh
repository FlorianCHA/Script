#!/bin/sh
# properties = {"threads": 4, "params": {}, "input": ["/homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta", "/homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa"], "jobid": 8, "log": [], "wildcards": ["AG0004"], "local": false, "resources": {}, "cluster": {}, "rule": "exonerate", "output": ["/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.gff3", "/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3"]}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.gff3 /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_AG0004.hints.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.kxf4aiy0 /homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta /homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.kxf4aiy0/8.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.kxf4aiy0/8.jobfailed"; exit 1)

