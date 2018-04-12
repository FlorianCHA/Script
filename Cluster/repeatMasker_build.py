#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package repeatMasker_build.py
# @author Florian Charriat

"""
	The repeatMasker_build script
	=============================
	:author: CHARRIAT Florian\n
	:contact: florian.charriat@inra.fr\n
	:date: 15/03/2018\n
	:version: 0.1\n

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
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFasta , recupId



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to create a bash file that can launch jobs on the cluster that uses repeatMasker on all the fasta files in a directory. ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display repeatMasker_build.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the fasta files that repeatMasker should use')
	filesreq.add_argument('-b', '--database',type = str, required=True, dest = 'database', help = 'Path of database file')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	database=  os.path.abspath(args.database)
	outDir= os.path.abspath(args.outdirPath)


########### Gestion directory #############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	verifFichier(database)
	bash = outDir+'script_bash'
	name_directory = [outDir, bash, outDir+'result_repeatMasker', outDir+'sge_output',outDir+'sge_error']
	for folder in name_directory: 
		createDir(folder)
		
############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in repeatMasker_build (Version " + version + ")     ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')
	
########## main script ######################
	nbScript = 0
	runJob = open(outDir+"run_repeatMAskerJob.sh","w")
	for files in os.listdir(directory):
		if isFasta(files) :
			nbScript += 1
			filesName = recupId(files)
			createDir(outDir +"result_repeatMasker/"+filesName)
			SCRIPT = open(bash+'/'+filesName+ "_repeatMasker.sh","w")
			SCRIPT.write('#$ -o '+outDir+'sge_output/'+filesName+'.out\n#$ -e '+outDir+'sge_error/'+filesName+'.err\nmodule load bioinfo/RepeatMasker/4.0.7;\n')
			SCRIPT.write("RepeatMasker -pa 4 -s -no_is -nolow "+directory+files+" -lib "+database+" -e ncbi -dir "+ outDir +"result_repeatMasker/"+filesName+";\n")
			SCRIPT.close()
			os.system("chmod 755 "+outDir+'script_bash/'+filesName+ "_repeatMasker.sh")
			runJob = open(outDir+"run_repeatMAskerJob.sh","a")
			runJob.write("qsub -N "+filesName+"_repeatmasker -V -q long.q -pe parallel_smp 4 "+ outDir+'script_bash/'+filesName+ "_repeatMasker.sh\n")
			runJob.close()
			print(filesName +'\t done')		
			

############## summary message #######################

	print(form('\n-------------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\tInput :')
	print('\t\t- Repertoire assemblage : '+directory[:-1])
	print('\t\t- Database : '+database)
	
	print('\n\tOutput :')
	print('\t\t- script bash : '+bash)
	print('\t\t- fichier a lancer : '+outDir+"run_repeatMAskerJob.sh")	
	print('\t\t- Resultat des Jobs : '+outDir +'result_repeatMasker')	
	
	print('\n Si vous souhaité lancer les '+str(nbScript)+' jobs veuillez taper la commande : ')
	print(form('\n\t\t\t\tbash '+outDir+'run_job_mapping.sh\n','green','bold'))
	print(form('-------------------------------------------------------------------------------------------------------------------------','red','bold'))

	

			
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))


