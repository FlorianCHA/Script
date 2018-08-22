#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package Recup_TO_Only.py
# @author Florian Charriat

"""
	The Recup_TO_Only script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 15/08/2018
	:version: 0.1

	Script description
	------------------

	This program is used to remove all alternative transcript of protein fasta file. 

	Example
	-------

	>>> Recup_TO_Only.py -d /homedir/user/work/data/ -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display Recup_TO_Only.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta files which must correct
						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""
######### Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''	This program is used to remove all alternative transcript of protein fasta file.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display ABYSS_launch version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the fasta files which must correct')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	output= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)
	output = verifDir(output)
	name_directory = [output]
	createDir(name_directory)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in Recup_To_Only (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')
	



	#directory = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/1_MGG/'
	#output = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/1_MGG/'

	for fasta in os.listdir(directory):
		print('%s in process'%fasta)
		dico = fasta2dict(directory+fasta)
		f = open(output+fasta,'w')
		for idSeq in sorted(dico.keys(), key=sort_human):
			if 'T0' in idSeq.split('_')[-1] :
				seqObj = dico[idSeq].seq
				record = SeqRecord(seqObj,id=idSeq,name=idSeq, description= dico[idSeq].description)
				SeqIO.write(record,f, "fasta")

		
