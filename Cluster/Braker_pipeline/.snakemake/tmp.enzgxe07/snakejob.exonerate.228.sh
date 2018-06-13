#!/bin/sh
# properties = {"rule": "exonerate", "cluster": {}, "log": [], "input": ["/homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BR29.fasta"], "params": {}, "wildcards": ["BR29"], "threads": 4, "resources": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_BR29.gff3", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_BR29.hints.gff3"], "jobid": 228, "local": false}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_BR29.gff3 /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_BR29.hints.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07 /homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/BR29.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07/228.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.enzgxe07/228.jobfailed"; exit 1)

