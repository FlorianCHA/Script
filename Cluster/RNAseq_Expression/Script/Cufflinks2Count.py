#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
import sys

pathfile = sys.argv[1]
RNAseq = sys.argv[2]
f = open(pathfile+'/transcripts.gtf','r')
lines = f.readlines()
f.close()

for line in lines :
	tabLine = line.split('\t')
	if tabLine[2] == 'transcript' :
		position = '%s:%s_%s'%(tabLine[0],tabLine[3],tabLine[4])
		tabScore = tabLine[8].split('"')
		name = tabScore[3]
		FPKM = tabScore[5]
		couverture = tabScore[-4]
		print('%s\t%s\t%s\t%s\t%s'%(name,position,FPKM,couverture,RNAseq))
