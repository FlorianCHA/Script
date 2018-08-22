#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package Verif_Annotation.py
# @author Florian Charriat

"""
	The Verif_Annotation script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to verif annotation (start or stop error or multi-stop)

	Example
	-------

	>>> Verif_Annotation.py -d /homedir/user/work/data -

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display Verif_Annotation.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains the result of annotation pipeline


"""


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to verif annotation (start or stop error or multi-stop)''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Verif_Annotation version number and exit')
	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains the result of annotation pipeline')


	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outDir= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash'
	name_directory = [outDir,outDir+'error_files', outDir+'out_files',bash,outDir+'result']
	createDir(name_directory)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in Verif_Annotation (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

	#directory = '/homedir/gladieux/work/magMax_project/2_Annotation/4_final-data/'

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
		
		

		
		
		
		
		
