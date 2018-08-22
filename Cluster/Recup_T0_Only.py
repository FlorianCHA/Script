#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-

########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens


directory = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/1_MGG/'
output = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/1_MGG/'

for fasta in os.listdir(directory):
	print('%s in process'%fasta)
	dico = fasta2dict(directory+fasta)
	f = open(output+fasta,'w')
	for idSeq in sorted(dico.keys(), key=sort_human):
		if 'T0' in idSeq.split('_')[-1] :
			seqObj = dico[idSeq].seq
			record = SeqRecord(seqObj,id=idSeq,name=idSeq, description= dico[idSeq].description)
			SeqIO.write(record,f, "fasta")

		
