#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The ABYSS_launch script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

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
						show this help 	message and exit
		- \-v, --version
						display ABYSS_launch.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta files which must be assembled
						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display ABYSS_launch version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the fasta files which must be assembled')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
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
	print("\t"+form("|",'yellow','bold')+form("        Welcome in ABYSS_launch (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################
	nbJob = 0
	nbGenome = 0
	for file in os.listdir(directory):
		if file.endswith('_R1.fastq.gz')==True :
			nbGenome += 1
			isolate=file.replace('_R1.fastq.gz','')
			print(form('\nLancement des jobs pour : '+isolate+'\n','green',['bold','underline']))
			for kmers in [20,30,40,50,60,70,80,90]:
				nameScript =bash+'/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh'
				SCRIPT=open(nameScript,'w')
				SCRIPT.write('#$ -o '+outDir+'out_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.out\n#$ -e '+outDir+'error_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.err\nmodule load bioinfo/abyss/1.9.0;\n')
				SCRIPT.write('mkdir -p '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write('cd '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write("/usr/local/bioinfo/abyss/1.9.0/bin/abyss-pe name="+isolate+"_"+str(kmers)+" k="+str(kmers)+" in='"+directory+file+" "+directory+file.replace('_R1','_R2')+"' -o abyss_assembly_"+isolate+"_"+str(kmers)+".fasta;\n")
				SCRIPT.close()
				os.system('qsub -l mem_free=50G -l h_vmem=60G -q normal.q '+nameScript)
				nbJob += 1


############## summary message #######################
	print(form('\n---------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\n\tInput : \n\t\t- '+ directory[:-1])
	print('\n\tOutput :')
	print('\t\t - script bash créé : ' +bash)
	print('\t\t - Résultats des assemblages : '+outDir+'result')
	print('\n\tAbyss_launch a lancés '+str(nbJob)+" jobs pour l'assemblage des "+str(nbGenome)+' souches\n')
	print(form('---------------------------------------------------------------------------------','red','bold'))

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




