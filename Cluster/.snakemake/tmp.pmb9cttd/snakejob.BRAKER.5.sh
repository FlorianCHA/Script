#!/bin/sh
# properties = {"output": ["/homedir/charriat/work/Annotation/test/2_Braker/AG0004/"], "log": [], "params": {}, "local": false, "jobid": 5, "cluster": {}, "threads": 4, "input": ["/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_AG0004.gff", "/homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta"], "wildcards": ["AG0004"], "resources": {}, "rule": "BRAKER"}
cd /gs7k1/home/charriat/Script/Cluster && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/2_Braker/AG0004/ --snakefile /gs7k1/home/charriat/Script/Cluster/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.pmb9cttd /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints.intron_AG0004.gff /homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules BRAKER  && touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.pmb9cttd/5.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/.snakemake/tmp.pmb9cttd/5.jobfailed"; exit 1)

