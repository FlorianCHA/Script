#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package Verif_Annotation.py
# @author Florian Charriat



########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human

directory = '/homedir/gladieux/work/magMax_project/2_Annotation/4_final-data/'

print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('Isolat','Nombre Met','NoStart','NoStop','MultiStop','NoStartAndStop','NoStartAndMultiStop','NoStopAndMultiStop','NoAll'))
for folder in os.listdir(directory) :
	Nostart = 0
	Nostop = 0
	MultiStop = 0
	NoStartAndStop = 0
	NoStartAndMultiStop = 0
	NoStopAndMultiStop = 0
	NoAll = 0
	NbMet = 0
	pathFile = directory+folder+'/'+folder+'_protein.faa'
	dico = fasta2dict(pathFile)
	for elt  in sorted(dico.keys(), key=sort_human):
		seq = str(dico[elt].seq)
		NbMet += seq.count('M')
		if seq[0] != 'M' :
			if seq[-1] != '*' :
				if seq.count('*') != 1:
					NoStartAndMultiStop += 1
				else :
					NoStartAndStop += 1
			
			if seq.count('*') != 1:
				NoAll += 1
			else : 
				Nostart += 1
			
			
		elif seq[-1] != '*' :
			if seq.count('*') != 1 :
				NoStopAndMultiStop += 1
			else :
				Nostop += 1			
		
		elif seq.count('*') != 1:
			MultiStop += 1
	print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(folder,NbMet,Nostart,Nostop,MultiStop,NoStartAndStop,NoStartAndMultiStop,NoStopAndMultiStop,NoAll))
		
		

		
		
		
		
		
