#!/bin/sh
# properties = {"local": false, "output": ["/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.gff3", "/homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.hints.gff3"], "params": {}, "rule": "exonerate", "cluster": {}, "input": ["/homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa", "/homedir/charriat/work/Annotation/test_input/CH1857_scaffold.fasta"], "wildcards": ["CH1857"], "threads": 4, "log": [], "resources": {}, "jobid": 7}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.gff3 /homedir/charriat/work/Annotation/test/1_hints/ProtHints/exonarate_CH1857.hints.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.l_nld1q5 /homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa /homedir/charriat/work/Annotation/test_input/CH1857_scaffold.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.l_nld1q5/7.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.l_nld1q5/7.jobfailed"; exit 1)

