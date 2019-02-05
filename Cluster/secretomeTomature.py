#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package secretomeTomature.py.py
# @author Florian CHARRIAT

"""
	The secretomeTomature.py script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to retrieve mature protein to secretome fasta file
	Example
	-------
	>>> secretomeTomature.py.py -d /homedir/secretome/ -o result/


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display secretomeTomature.py.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory which contains all secretome fasta file
		- \-o <path/to/directory/output>, --output <path/to/directory/output>
						path of the output directory

"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human

## Python modules
import argparse


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This programme is used to retrieve mature protein to secretome fasta file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>",type = str, required=True, dest = 'dirPath', help = 'path of directory which contains all secretome fasta file')
	filesreq.add_argument('-o', '--output', metavar="<path/to/directory/clumpak>",type = str, required=True, dest = 'output', help = 'path of the output directory')


######### Recuperation arguments ###########
	args = parser.parse_args()
	dir = os.path.abspath(args.dirPath)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	dir = verifDir(dir,True)
	output = verifDir(output)
	createDir([output])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in vcfSearch (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	with open(f'{output}Launch_secretome_to_mature_script.sh','w') as f :
		for elt in os.listdir(dir):
			filePath = f'{dir}{elt}'
			outPath = f'{output}{elt.replace(".fasta","")}_mature.fasta'
			with open(f'{output}{elt.replace(".fasta","")}.sh','w') as f1:
				f1.write(f'signalp -m {outPath} -u 0.01 -U 0.01 {filePath}')
			f.write(f'qsub -q normal.q -V -cwd {output}{elt.replace(".fasta","")}.sh\n')


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))