#!/bin/env python
# -*- coding: utf-8 -*-
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
## Python modules
import argparse, os, sys

#Import MODULES_SEB
from module_Flo import verifDir, createDir



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
				SCRIPT.write('# !/bin/bash\n#$ -o '+outDir+'out_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.out\n#$ -e '+outDir+'error_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.err\n#$ -N '+isolate+'_'+str(kmers)+'\n#$ -l h_vmem=36G\n#$ -l mem=32G\n')
				SCRIPT.write('mkdir -p '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write('cd '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write("/usr/local/bioinfo/src/abyss/abyss-2.0.0/bin/abyss-pe name="+isolate+"_"+str(kmers)+" k="+str(kmers)+" in='"+directory+file+" "+directory+file.replace('_R1','_R2')+"' -o abyss_assembly_"+isolate+"_"+str(kmers)+".fasta;\n")
				SCRIPT.close()
				os.system('chmod 755 '+outDir+'script_bash/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh')
				os.system('qsub "'+outDir+'script_bash/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh"')

