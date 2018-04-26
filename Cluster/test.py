#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package test.py
# @author Florian Charriat




from module_Flo import parseGFF

files = 'CH1857'
out_file = '/homedir/charriat/work/Annotation/test_braker/Braker/%s/braker/magnaporthe_oryzae/test.gff3'%files
GFFfile = parseGFF('/homedir/charriat/work/Annotation/test_braker/Braker/%s/braker/magnaporthe_oryzae/augustus.gff3'%files)





with open(out_file, "w") as out_handle:
	for record in GFFfile.parseGFF3():
		test = []
		for elt in record.attributes :

			record.attributes[elt] = record.attributes[elt].replace('g',files+'_')
			test.append('%s=%s'%(elt,record.attributes[elt]))
		if len(test) == 1 :
			test = test[0]
		if len(test) == 2 :
			test = '%s;%s'%(test[0],test[1])
		out_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(record.seqid,record.source,record.type,record.start,record.end,record.score,record.strand,record.phase,test))
