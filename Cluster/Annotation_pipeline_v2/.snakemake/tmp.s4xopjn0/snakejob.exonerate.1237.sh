#!/bin/sh
# properties = {"threads": 2, "local": false, "input": ["/work/gladieux/magMax_project/2_New_Annotation/0_rawdata/CH1876_scaffold.fasta", "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/SupplementaryFile/OG_filter.fasta"], "log": [], "cluster": {}, "rule": "exonerate", "wildcards": ["CH1876_scaffold"], "jobid": 1237, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_CH1876_scaffold.hints.gff3", "/work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_CH1876_scaffold.gff3"], "params": {"l_mem_free": "4G"}, "resources": {}}
cd /work/gladieux/Script/Cluster/Annotation_pipeline_v2 && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_CH1876_scaffold.hints.gff3 /work/gladieux/magMax_project/2_New_Annotation/1_tmp/1_hints/ProtHints/exonerate_CH1876_scaffold.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline_v2/BRAKER_pipeline2.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0 /work/gladieux/magMax_project/2_New_Annotation/0_rawdata/CH1876_scaffold.fasta /work/gladieux/Script/Cluster/Annotation_pipeline_v2/SupplementaryFile/OG_filter.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/1237.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.s4xopjn0/1237.jobfailed"; exit 1)

