#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package selection_TMHMM.py
# @author Florian Charriat

"""
	The selection_TMHMM script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/04/2018
	:version: 0.1

	Script description
	------------------

	This program is used to selected the proteins having no transmembrane domain or only one in the 60 first amino acid. The program take in input the output of TMHMM Tools. This program is used by secretome_Pipeline.
	
	Example
	-------

	>>> selection_TMHMM.py  -f /homedir/user/work/data/TMHMM_result.txt -o /homedir/user/work/result/


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display selection_TMHMM.py version number and exit
						
	Input mandatory infos for running:
								
		- \-t <path/to/TMHMM/output/file>, --file <path/to/TMHMM/output/file>
						path of the TMHMM output file
		- \-f <path/to/fasta/file>, --file <path/to/fasta/file>
						Path of the fasta file which TMHMM has been proceed
		- \-o <path/to/output/file>, --outdirPath <path/to/output/file>
						path and name of the output file

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
	parser = argparse.ArgumentParser(prog=__file__, description='''
	This program is used to selected the proteins having no transmembrane domain or only one in the 60 first amino acid. The TMHMM Tools is used. This program is used by secretome_Pipeline.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display selection_TMHMM version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-t', '--TMHMM',type = str, default = 'None', dest = 'TMHMM', help = 'Path of the TMHMM output file')
	filesreq.add_argument('-f', '--fasta',type = str,  required=True, dest = 'fasta', help = 'Path of the fasta file which TMHMM has been proceed')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdir', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	outDir= os.path.abspath(args.outdir)
	TMHMM = os.path.abspath(args.TMHMM)
	fasta = os.path.abspath(args.fasta)
	Id = recupId(fasta.split('/')[-1])
	
########### Gestion directory ##############

	name_directory = [outDir]
	createDir(name_directory)
	outDir = verifDir(outDir)

	


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in selection_TMHMM   (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')
	
	
	
	
##################### Main ###########################

	inputFile = open(TMHMM,'r')
	lines = inputFile.readlines()
	nbFaux = 0
	nbSecretome = 0
	nbMoyen = 0
	listeID = []
	text = ''

	for line in lines :
		lineSplit = line.split('\t')
		if lineSplit[4] == 'PredHel=0' :
			nbSecretome += 1
			text = text + line
			listeID.append(lineSplit[0])
		elif lineSplit[4] == 'PredHel=1' :
			FirstTH_s = int(lineSplit[5].split('=')[1].split('-')[0].replace('o','').replace('i',''))
			FirstTH_e = int(lineSplit[5].split('=')[1].split('-')[1].replace('o','').replace('i',''))
			if  FirstTH_s < 40 and FirstTH_e < 40 :
				nbSecretome += 1
				text = text + line
				listeID.append(lineSplit[0])
			elif FirstTH_s < 40 and FirstTH_e > 40 :
				nbMoyen +=1 
		else :
			nbFaux += 1
	
	f = open('%s%s_selectTMHMM.txt'%(outDir,Id),'w')
	f.write('%s\n%s%s%s\n%s\n\nTotal processed sequence : %s\nNumber of protein with a single TH in the first 60 aa : %s\nnNumber of protein with a single TH not only the first 60 aa : %s\nEliminer : %s\n\n'%('#'*38,'#'*10,' Protein selected ','#'*10,'#'*38,nbFaux+nbMoyen+nbSecretome,nbSecretome,nbMoyen,nbFaux))
	f.write(text)
	f.close()
	
	dico_fasta = fasta2dict(fasta)
	f = open('%s%s_secreted_2.fasta'%(outDir,Id),'w')
	
	for idSeq in sorted(listeID, key=sort_human):
			seqObj = dico_fasta[idSeq].seq
			record = SeqRecord(seqObj,id=idSeq,name=idSeq, description=dico_fasta[idSeq].description)
			SeqIO.write(record,f, "fasta")

############## summary message #######################


	print(form('\n-----------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\n\tInput : \n\t\t- '+ fasta+'\n\t\t-'+TMHMM)
	print('\n\tOutput :')
	print('\t\t - Résultat des prédictions des secretomes : '+outDir)
	print('%s\n\t%s%s%s\n\t%s\n\n\tTotal processed sequence : %s\n\tNumber of protein with a single TH in the first 60 aa : %s\n\tNumber of protein with a single TH not only the first 60 aa : %s\n\tEliminer : %s\n'%('#'*38,'#'*10,' Protein selected ','#'*10,'#'*38,nbFaux+nbMoyen+nbSecretome,nbSecretome,nbMoyen,nbFaux))
	print(form('----------------------------------------------------------------------------------------------------------------------','red','bold'))

		

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))





