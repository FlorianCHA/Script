#!/usr/bin/python2.7
#-*- coding: utf-8 -*-
# @package repeatMasker_build.py
# @author Florian Charriat

"""
	The repeatMasker_build script
	=============================
	author: CHARRIAT Florian\n
	contact: florian.charriat@inra.fr\n
	date: 15/03/2018\n
	version: 0.1\n

	Script description
	------------------

	This program is used to create a bash file that can launch jobs on the cluster that uses repeatMasker on all the fasta files in a directory.

	Example
	-------

	>>> repeatMasker_build.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display repeatMasker_build.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta files that repeatMasker should use
		- \-b <path/to/database>, --database <path/to/database>
						path of database file
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to create a bash file that can launch jobs on the cluster that uses repeatMasker on all the fasta files in a directory. ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the fasta files that repeatMasker should use')
	filesreq.add_argument('-b', '--database',type = str, required=True, dest = 'database', help = 'path of database file', default = '/gs7k1/projects/BGPI/becphy/pangenome2017/fungi_refTE70-15.fasta')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	database=  args.database
	outDir= args.outdirPath


########### Gestion directory #############
	directory = verifDir(directory)
	database = database
	outDir = verifDir(outDir)

	name_directory = [outDir, outDir+'script_bash', outDir+'result_repeatMasker', outDir+'sge_output',outDir+'sge_error']
	for folder in name_directory: 
		createDir(folder)

########## main script ######################
	runJob = open("run_repeatMAskerJob.sh","w")
	for files in os.listdir(directory):
		if files.endswith('.fasta')==True or files.endswith('.fa')==True :
			filesName = files.replace('.fa','')
			filesName = filesName.replace('.fasta','')
			createDir(outDir +"result_repeatMasker/"+filesName)
			SCRIPT = open(outDir+'script_bash/'+filesName+ "_repeatMasker.sh","w")
			SCRIPT.write('#$ -o '+outDir+'sge_output/'+filesName+'.out\n#$ -e '+outDir+'sge_error/'+filesName+'.err\nmodule load bioinfo/RepeatMasker/4.0.7;\n')
			SCRIPT.write("RepeatMasker -pa 4 -s -no_is -nolow "+directory+files+" -lib "+database+" -e ncbi -dir "+ outDir +"result_repeatMasker/"+filesName+";\n")
			SCRIPT.close()
			os.system("chmod 755 "+outDir+'script_bash/'+filesName+ "_repeatMasker.sh")
			runJob = open("run_repeatMAskerJob.sh","a")
			runJob.write("qsub -N "+filesName+"_repeatmasker -V -q long.q -pe parallel_smp 4 "+ outDir+'script_bash/'+filesName+ "_repeatMasker.sh\n")
			runJob.close()
			print(files +'\t done')		

