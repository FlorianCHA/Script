#!/bin/sh
# properties = {"threads": 2, "input": ["/homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta", "/homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_AG0004.gff"], "log": [], "output": ["/homedir/charriat/work/Annotation/test/3_augustus/AG0004_augustus.gff3"], "wildcards": ["AG0004"], "params": {}, "resources": {}, "local": false, "cluster": {}, "jobid": 5, "rule": "augustus"}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /homedir/charriat/work/Annotation/test/3_augustus/AG0004_augustus.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.p9twb9vc /homedir/charriat/work/Annotation/test_input/AG0004_scaffold.fasta /homedir/charriat/work/Annotation/test/1_hints/MergeHints/RNAseq_protein.hints_AG0004.gff --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules augustus  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.p9twb9vc/5.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.p9twb9vc/5.jobfailed"; exit 1)

