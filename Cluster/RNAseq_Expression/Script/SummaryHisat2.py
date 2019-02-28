#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
import sys

pathfile = sys.argv[1]
RNAseq = sys.argv[2]
f = open(pathfile,'r')
lines = f.readlines()
f.close()

for line in lines :
	if 'reads; of these:' in line :
		nbAll = line.split()[0]
	if '0 times' in line :
		noAln = line.strip().split()[0]
	if 'exactly 1 time' in line :
		Aln1 = line.strip().split()[0]
	print('%s\t%s\t%s\t%s\t%s' % (nbAll, noAln, Aln1,RNAseq))
