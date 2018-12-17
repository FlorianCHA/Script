#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package seachHmm.py
# @author Florian Charriat

"""
	The seachHmm script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to search sequence with a Hmm profil and create a fasta file.

	Example
	-------

	>>> seachHmm.py -f /homedir/user/work/fasta.file -p /homedir/user/work/profil -o /homedir/user/work/result.fasta

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display seachHmm.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta that contains all the sequence to use.
						
		- \-p <path/to/HMM/profil/file>, --profil <path/to/HMM/profil/file>
						path of the Hmm profil which be used

		- \-e <path/to/HMM/profil/file>, --evalue <path/to/HMM/profil/file>
						Expect E-value for saving hits (default = 10). For exemple you can put 1e-4 for 10^-4

		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file (fasta file)

"""


########## Module ###############
## Python modules
import argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form,verifFichier,fasta2dict



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to seach sequence with a Hmm profil and create a fasta file.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display seachHmm version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of fasta that contains all the sequence to use.')
	filesreq.add_argument('-e', '--evalue',type = float, required=False,default = 10, dest = 'evalue', help = 'Expect E-value for saving hits (default = 10). For exemple you can put 1e-4 for 10^-4.\n ')
	filesreq.add_argument('-p', '--profil',type = str, required=True, dest = 'profil', help = 'Path of the Hmm profil which be used')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	profil = os.path.abspath(args.profil)
	output = os.path.abspath( args.output)
	evalue = args.evalue

########### Gestion directory ##############
	verifFichier(fasta)
	verifFichier(profil)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("            Welcome in seachHmm (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

########### Main #####################################
		
	result_output = output.replace('.fasta','_result.txt')
	print(form('Recherche de séquence avec Hmmsearch\n\n\t- Profil utilisé : %s\n\t- Fasta utilisé : %s\n\t- Evalue : %s\n'%(profil,fasta,evalue),'white','bold'))
	dico_fasta = fasta2dict(fasta)
	os.system('hmmsearch  --max --nonull2 -E %s %s %s > %s'%(evalue,profil,fasta,result_output))

	print(form("---------------------------------------------------------------------------------------------------------------------------------",'yellow','bold')+'\n')

	print(form('Création du fichier fasta\n\n\t- fichier résultat : %s\n\t- output : %s\n'%(result_output,output),'white','bold'))
	file_result = open(result_output,'r')
	lines = file_result.readlines()
	file_result.close()
	start = False
	recup = False
	liste = []
	output_file = open(output,'w')
	for line in lines :
		line = line.strip()
		if line[0:7] == 'E-value' :
			start = True 
		elif '------' in line and start == True :
			recup = True 
		elif recup == True and line == ''  :
			break
		elif recup == True :
			lineSplit = line.split()
			e_full = lineSplit[0]
			e_domain = lineSplit[3]
			name = lineSplit[8]
			seq = dico_fasta[name].seq
			length = len(str(dico_fasta[name].seq))
			record = SeqRecord(seq,id=name,name=name, description=' | evalue full sequence = %s | evalue best domain = %s length = %s'%(e_full,e_domain,length)) 
			SeqIO.write(record,output_file, "fasta")

	output_file.close()



############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------\n",'yellow','bold'))



















