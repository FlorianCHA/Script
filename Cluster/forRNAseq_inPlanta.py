#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The forRNAseq_inPlanta script
	=============================

	author: CHARRIAT Florian
	contact: florian.charriat@inra.fr
	date: 21/03/2018
	version: 0.1

	Script description
	------------------

	This program is used to align multiple RNA-seq on a host genome for kept only RNA-seq of interest

	Example
	-------

	>>> forRNAseq_inPlanta .py -d /homedir/user/work/dataRNAseq -o /homedir/user/work/result -r /homedir/user/work/genome.fasta -c /homedir/user/work/configfile.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display ABYSS_launch.py version number and exit
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
import argparse
import os


####### FUNCTION ################
def createDir(directory):
	''' Permet de vérifier si un dossier existe, si ce n'est pas le cas, 
	le dossier sera crée
	'''
	if not os.path.exists(directory):
		 	os.makedirs(directory)
	return

def verifDir(directory):
	'''
	Permet de mettre en forme le chemin du dossier pour être utilisé dans un script, 
	la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute
	'''
	if directory.endswith('/') == False :
		directory = directory +'/'

	return directory


if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to align multiple RNA-seq on a genome and merge the different alignments into a file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the RNA-seq')
	filesreq.add_argument('-r', '--refDir',type = str, required=True, dest = 'refDir', help = 'path of directory that contains  the host genome ')
	filesreq.add_argument('-c', '--configFile',type = str, required=True, dest = 'configFile', help = 'Ppath of the config file for toogle')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	genome = args.refDir
	config = args.configFile
	outDir= args.outDirPath


########### Gestion directory ##############
	directory = verifDir(directory)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash/'
	name_directory = [outDir, bash,outDir+'error/',outDir+'output/']
	for folder in name_directory: 
		createDir(folder)

########## Main ############################
	IDgenome = genome.split('/')
	IDgenome = IDgenome[-1].replace('.fasta','')
	genomeOutDir = outDir+IDgenome
	nameFile = bash+IDgenome+'_mapping.sh'
	files = open(nameFile,'w')
	#Permet de géré les sortie de sge 
	files.write('#$ -o '+outDir+'output/'+IDgenome+'.out\n#$ -e '+outDir+'error/'+IDgenome+'.err\n')
	# Permet de charger puis lancer Toogle pour un alignement
	files.write('module load bioinfo/braker/1.9 bioinfo/exonerate/2.4.7 bioinfo/TOGGLE/dev bioinfo/R/3.2.2\n')
	files.write('toggleGenerator.pl -d '+directory+' -r '+genome+' -c '+config+' -o '+genomeOutDir+' -nocheck;\n')
	files.close()	
	os.system("qsub -N "+IDgenome+"_mapping -V -q normal.q '"+nameFile+"'\n")









