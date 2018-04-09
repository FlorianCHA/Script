#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The forRNAseq_inPlanta script
	=============================

	:author: CHARRIAT Florian
	:contact: florian.charriat@inra.fr
	:date: 21/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to align multiple RNA-seq on a host genome for eliminate inPlanta RNA

	Example
	-------

	>>> forRNAseq_inPlanta .py -d /homedir/user/work/dataRNAseq -o /homedir/user/work/result -r /homedir/user/work/genome.fasta -c /homedir/user/work/configfile.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display forRNAseq_inPlanta.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the RNA-seq in planta
		- \-r <path/to/file/reference/>, --refDir<path/to/reference/directory>
						path of directory that contains the host genome
		- \-c <path/to/configFile>, --configFile <path/to/>configfile>
						path of the config file for toogle
		- \-o <path/to/output/directory>, --outDirPath <path/to/output/directory>
						Path of the output directory


"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFastq , recupId



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to align multiple RNA-seq on a host genome for eliminate inPlanta RNA''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the RNA-seq')
	filesreq.add_argument('-r', '--refDir',type = str, required=True, dest = 'refDir', help = 'path of directory that contains  the host genome ')
	filesreq.add_argument('-c', '--configFile',type = str, required=True, dest = 'configFile', help = 'Ppath of the config file for toogle')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	genome = os.path.abspath(args.refDir)
	config = args.configFile
	outDir= os.path.abspath(args.outDirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	verifFichier(config)
	bash = outDir+'script_bash/'
	name_directory = [outDir, bash,outDir+'error/',outDir+'output/']
	for folder in name_directory: 
		createDir(folder)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("     Welcome in forRNAseq_inPlanta (Version " + version + ")       ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

########## Main ############################
	nbRNAseq = 0
	for data in os.listdir(directory) :
		if isFastq(data):
			nbRNAseq += 1
	IDgenome = genome.split('/')
	IDgenome = recupId(IDgenome[-1])
	genomeOutDir = outDir+IDgenome
	nameFile = bash+IDgenome+'_mapping.sh'
	files = open(nameFile,'w')
	#Permet de géré les sortie de sge 
	files.write('#$ -o '+outDir+'output/'+IDgenome+'.out\n#$ -e '+outDir+'error/'+IDgenome+'.err\n')
	# Permet de charger puis lancer Toogle pour un alignement
	files.write('module load bioinfo/TOGGLE/0.3.6\n')
	files.write('rm -rf '+genomeOutDir+'\n') # permet de suprimer le fichier output de TOGGLE pour éviter tous problèmes
	files.write('toggleGenerator.pl -d '+directory[:-1]+' -r '+genome+' -c '+config+' -o '+genomeOutDir+' -nocheck;\n')
	files.close()	
	os.system("qsub -N "+IDgenome+"_mapping -V -q normal.q '"+nameFile+"'\n")


############## summary message #######################
	print(form('\n---------------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	
	print('\tInput :')
	print('\t\t- Donnée RNAseq : '+directory[:-1])
	print('\t\t- Repertoire génome : '+genome)
	print('\t\t- Fichier config : '+config)
	
	print('\n\tOutput :')
	print('\t\t- script bash : '+bash)
	print('\t\t- Resultat du Job : '+ genomeOutDir)	
	
	
	
	
	
	print('\n forRNAseq_inPlanta a lancés le job de mapping des '+str(nbRNAseq)+' données de RNAseq InPlanta sur le genome de reférence ('+IDgenome+')\n')

	print(form('----------------------------------------------------------------------------------------------------------------------------','red','bold'))

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))







