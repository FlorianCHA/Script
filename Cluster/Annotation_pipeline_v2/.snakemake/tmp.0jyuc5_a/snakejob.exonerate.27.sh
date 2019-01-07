#!/bin/sh
# properties = {"input": ["/work/gladieux/magMax_project/pacbio_canu_arrow/PH0014-arrow2.fasta", "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/SupplementaryFile/OG.fasta"], "resources": {}, "rule": "exonerate", "cluster": {}, "wildcards": ["PH0014-arrow2"], "jobid": 27, "params": {"l_mem_free": "4G"}, "local": false, "output": ["/work/gladieux/magMax_project/2_New_Annotation/Canu-Annotation/1_hints/ProtHints/exonerate_PH0014-arrow2.gff3", "/work/gladieux/magMax_project/2_New_Annotation/Canu-Annotation/1_hints/ProtHints/exonerate_PH0014-arrow2.hints.gff3"], "threads": 2, "log": []}
cd /work/gladieux/Script/Cluster/Annotation_pipeline_v2 && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/Canu-Annotation/1_hints/ProtHints/exonerate_PH0014-arrow2.gff3 /work/gladieux/magMax_project/2_New_Annotation/Canu-Annotation/1_hints/ProtHints/exonerate_PH0014-arrow2.hints.gff3 --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline_v2/BRAKER_pipeline2.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.0jyuc5_a /work/gladieux/magMax_project/pacbio_canu_arrow/PH0014-arrow2.fasta /work/gladieux/Script/Cluster/Annotation_pipeline_v2/SupplementaryFile/OG.fasta --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules exonerate  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.0jyuc5_a/27.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline_v2/.snakemake/tmp.0jyuc5_a/27.jobfailed"; exit 1)

