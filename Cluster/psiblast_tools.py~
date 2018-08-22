#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package psiblast_tools.py
# @author Florian Charriat

"""
	The psiblast_tools script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 4/06/2018
	:version: 0.1

	Script description
	------------------

	This program is used to align sequence with psiblast tools and give in supplementary a fasta file.

	Example
	-------

	>>> psiblast_tools.py -d /work/database/databaseName -f /work/database/query.fasta -o /work/result.fasta -i 10 -opt alignement -n exemple

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display psiblast_tools.py version number and exit
						
	Input  infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta file which must be used for alignement (query)
		- \-d <path/to/database>, --database <path/to/database/>
						path of database which we must be used for alignement
		- \-i <int>, --num_iteration <int>
						number of iteration for psiblast
		 \-opt <str>, --option <str>
						option for fasta file :
							\- A : aligned sequences
							\- G : aligned sequences with Gap
							\- C : complete sequence
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output fasta directory
		- \-n <str>, --name <str>
						name of output file (default = 'psiblast')
		- \-eval <str>, --evalue <float>
						Expect E-value for saving hits (default = 10). For exemple you can put 1e-4 for 10^-4.

"""


########## Module ###############
## Python modules
import argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from argparse import RawTextHelpFormatter
#Import module_Flo
from module_Flo import verifDir, createDir , form, fasta2dict



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to align sequence with psiblast tools and give in supplementary a fasta file. ''', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display psiblast_tools version number and exit')


	filesreq = parser.add_argument_group('Input infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'query', help = 'Path of fasta file which must be used for alignement (query)')
	filesreq.add_argument('-d', '--database',type = str, required=True, dest = 'database', help = 'Path of database which we must be used for alignement')
	filesreq.add_argument('-i', '--num_iteration',type = int, required=True, dest = 'interation', help = 'Number of iteration for psiblast')
	filesreq.add_argument('-opt', '--option',type = str, choices=['A','G','C'] ,required=True, dest = 'option', help = 'Option for fasta file : \n- A : aligned sequences \n- G : aligned sequences with Gap \n- C : complete sequence')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDir', help = 'Path of the output directory')
	filesreq.add_argument('-n', '--name',type = str, required=False,default = 'psiblast', dest = 'name', help = '\nName of output file (default = psiblast)\n ')
	filesreq.add_argument('-eval', '--evalue',type = float, required=False,default = 10, dest = 'evalue', help = 'Expect E-value for saving hits (default = 10). For exemple you can put 1e-4 for 10^-4.\n ')


	
######### Recuperation arguments ###########
	args = parser.parse_args()
	query = os.path.abspath(args.query)
	db = os.path.abspath(args.database)
	nbIteration = args.interation
	opt = args.option
	outDir= os.path.abspath(args.outDir)
	fileName = args.name
	evalue = args.evalue

########### Gestion directory ##############
	dbName = db.split('/')[-1]
	dbDirectory = db.replace(dbName,'')
	dbDirectory = verifDir(dbDirectory,True)
	outDir = verifDir(outDir)
	name_directory = [outDir]
	createDir(name_directory)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in psiblast_tools (Version " + version + ")        ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################.
	alignementFile = outDir+fileName+'_alignement'
	
############# Message Launch psiblast ######################
	print(form('Lancement de psiblast','green','bold'))
	print(form("-"*len('Lancement de psiblast'),'yellow','bold'))
	print(form('\nParamètres utilisés :',type='bold'))
	print('\n\t- Database : %s'%db)
	print('\t- Query : %s'%query)
	print('\t- Number of iteration : %s'%nbIteration)
	print('\t- Evalue min : %s'%evalue)
	print('\t- Output : %s\n\n'%alignementFile)
############### Fin message ################################

	os.system('psiblast -query %s -db %s -num_iterations %s -outfmt "6 qseqid sseqid pident evalue sseq" -evalue %s -out %s'%(query,db,nbIteration,evalue,alignementFile))
	f = open(alignementFile,'r')
	lines = f.readlines()
	f.close()
	
	############# Message creation fasta ######################
	print(form('Création fichier fasta','green','bold'))
	print(form("-"*len('Création fichier fasta'),'yellow','bold'))
	print(form('\nFasta type :',type='bold'))
	if opt == 'A':
		print('\n\t- Aligned sequences')
	if opt == 'G':
		print('\n\t- Aligned sequences with Gap')
	if opt == 'C':
		print('\n\t- Complete sequences')
	print(form('\nFichiers créés :',type='bold'))
	################ Fin message ##############################
		
	dico = fasta2dict(db+'.fasta')
	idQuery = 'None'
	for line in lines :
		if idQuery == 'None' :
			print('\n\t- %s_%s.fasta'%(fileName,line.split('\t')[0]))
			fastaFile = open('%s_%s.fasta'%(outDir+fileName,line.split('\t')[0]),"w")
		elif line.split('\t')[0] != idQuery :
			fastaFile.close()
			print('\t- %s_%s.fasta'%(fileName,line.split('\t')[0]))
			fastaFile = open('%s_%s.fasta'%(outDir+fileName,line.split('\t')[0]),"w")
		idQuery = line.split('\t')[0]
		idMatch = line.split('\t')[1]
		IndentityPer = line.split('\t')[2]
		Evalue = line.split('\t')[3]
		seqAlign = line.split('\t')[4]
		if opt == 'A':
			seq = seqAlign.replace('-','')[:-1]
			length = len(seq)
		if opt == 'G' :
			seq = seqAlign[:-1]
			length = len(seq)
		if opt == 'C':
			seq = str(dico[idMatch].seq)
			length = len(seq)
		record = SeqRecord(Seq(seq),id=idMatch,name=idMatch, description=' | evalue = %s | Percent_Identity = %s length = %s'%(Evalue,IndentityPer,length)) 
		SeqIO.write(record,fastaFile, "fasta")
		
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------\n",'yellow','bold'))
		
	












