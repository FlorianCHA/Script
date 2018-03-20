#!/usr/bin/python2.7
#-*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The ABYSS_launch script
	=======================
	author: Charriat Florian\n
	contact: florian.charriat@inra.fr\n
	date: 9/03/2018\n
	version: 0.1

	Script description
	------------------

	This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90)

	Example
	-------

	>>> ABYSS_launch.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display ABYSS_launch.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta files which must be assembled
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						Path of the output directory

"""


########## Module ###############

import argparse
import os

######## Fonction ###############
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the fasta files which must be assembled')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	outDir= args.outdirPath


########### Gestion directory ##############
	directory = verifDir(directory)
	outDir = verifDir(outDir)

	name_directory = [outDir,outDir+'error_files', outDir+'out_files',outDir+'script_bash',outDir+'result']
	for folder in name_directory: 
		createDir(folder)

############# Main #########################
	for file in os.listdir(directory):
		if file.endswith('_R1.fastq.gz')==True :
			isolate=file.replace('_R1.fastq.gz','')
			print(isolate)
			for kmers in [20,30,40,50,60,70,80,90]:
				SCRIPT=open(outDir+'script_bash/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh','w')
				SCRIPT.write('#$ -o '+outDir+'out_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.out\n#$ -e '+outDir+'error_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.err\nmodule load bioinfo/abyss/1.9.0;\n')
				SCRIPT.write('mkdir -p '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write('cd '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write("/usr/local/bioinfo/abyss/1.9.0/bin/abyss-pe name="+isolate+"_"+str(kmers)+" k="+str(kmers)+" in='"+directory+file+" "+directory+file.replace('_R1','_R2')+"' -o abyss_assembly_"+isolate+"_"+str(kmers)+".fasta;\n")
				SCRIPT.close()
				os.system('qsub -l mem_free=50G -l h_vmem=60G -q normal.q '+outDir+'script_bash/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh')

