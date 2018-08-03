#!/bin/sh
# properties = {"rule": "report", "input": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/data_report/Annotation_stat.csv", "/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/data_report/Assembly_quality.csv"], "threads": 1, "params": {}, "log": [], "output": ["/work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/report.html"], "wildcards": [], "resources": {}, "jobid": 38, "cluster": {}, "local": false}
cd /work/gladieux/Script/Cluster/Braker_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/report.html --snakefile /work/gladieux/Script/Cluster/Braker_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.xtuu5h9t /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/data_report/Annotation_stat.csv /work/gladieux/magMax_project/2_Annotation/5_soucheToulouse/6_report/data_report/Assembly_quality.csv --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules report  && touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.xtuu5h9t/38.jobfinished" || (touch "/work/gladieux/Script/Cluster/Braker_pipeline/.snakemake/tmp.xtuu5h9t/38.jobfailed"; exit 1)

