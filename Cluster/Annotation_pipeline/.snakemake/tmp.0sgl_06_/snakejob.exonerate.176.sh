#!/bin/sh
# properties = {"resources": {}, "local": false, "input": ["/work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/VT0030.fasta", "/homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa"], "wildcards": ["VT0030"], "params": {}, "cluster": {}, "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_VT0030.gff3", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_VT0030.hints.gff3"], "threads": 4, "rule": "exonerate", "log": [], "jobid": 176}
cd /gs7k1/home/charriat/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_VT0030.gff3 /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/1_hints/ProtHints/exonarate_VT0030.hints.gff3 --snakefile /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_ /work/gladieux/magMax_project/2_Annotation/0_rawdata/assembly_Toulouse/VT0030.fasta /homedir/charriat/BGPI/becphy/pangenome2017/test70-15/exonerate/70-15_annotated_protein.fa --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_/176.jobfinished" || (touch "/gs7k1/home/charriat/Script/Cluster/Braker_pipeline/.snakemake/tmp.0sgl_06_/176.jobfailed"; exit 1)

