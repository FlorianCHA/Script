#!/bin/sh
# properties = {"rule": "bamAlignement", "local": false, "resources": {}, "output": ["/work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1898_scaffold/", "/work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1898_scaffold/finalResults/bamList"], "cluster": {}, "params": {"l_mem_free": "4G"}, "jobid": 1554, "input": ["/homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/CH1898_scaffold.fasta", "/work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq", "SupplementaryFile/tophatMapping.config.txt"], "wildcards": ["/CH1898_scaffold"], "log": ["log/bamAlignement_{wildcards.smp}.out"], "threads": 1}
cd /work/gladieux/Script/Cluster/Annotation_pipeline && \
/gs7k1/binaries/snakemake/3.13.3/.venv/bin/python3 -m snakemake /work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1898_scaffold/ /work/gladieux/magMax_project/2_New_Annotation/1_tmp0_bamAlignement//CH1898_scaffold/finalResults/bamList --snakefile /work/gladieux/Script/Cluster/Annotation_pipeline/BRAKER_pipeline.snake \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg /homedir/gladieux/work/magMax_project/2_New_Annotation/0_rawdata/CH1898_scaffold.fasta /work/gladieux/magMax_project/2_Annotation/0_rawdata/rnaseq SupplementaryFile/tophatMapping.config.txt --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules bamAlignement  && touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg/1554.jobfinished" || (touch "/work/gladieux/Script/Cluster/Annotation_pipeline/.snakemake/tmp._vd2jvqg/1554.jobfailed"; exit 1)

