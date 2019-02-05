#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package parse_Busco.py
# @author Florian CHARRIAT

"""
	The parse_Busco script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme ....
	Example
	-------
	>>> parse_Busco.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display parse_Busco.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of result structure
		- \-c <path/to/directory/clumpak>, --clumpak <path/to/directory/clumpak>
						path of clumpak directory
		- \-l <filename>, --label <filename>
						File with LABEL, first column name, second top label
						info
	Input infos for running with default values:
		- \-dp <filename>, --drawparams <filename>
						Check your own drawparams file
		- \-co <filename>, --color <filename>
						File with colors (default 15 color max)
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme ....''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', type = str, required=True, dest = 'dir', help = 'path of directory which contain all Busco directory')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of output file')


######### Recuperation arguments ###########
	args = parser.parse_args()
	dir = os.path.abspath(args.dir)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	dir = verifDir(dir,True)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in vcfSearch (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	with open(output,'w') as output_file :
		output_file.write('Isolat\tComplete (%)\tFragmented (%)\tMissing (%)\n')
		for elt in os.listdir(dir) :
			if elt[0:4] =='run_' :
				isolat = elt.replace('run_','').replace('.busco','')
				with open(f'{dir}{elt}/short_summary_{isolat}.busco.txt','r') as summary_file:
					for line in summary_file:
						line = line.strip()
						if line[0:2] == 'C:' and 'n:290' in line :
							complete = line.split('%')[0].replace('C:','').replace('-','')
							imcompete = line.split('F:')[1].split('%')[0].replace('-','')
							miss = line.split('M:')[1].split('%')[0].replace('-','')
							output_file.write(f'{isolat}\t{complete}\t{imcompete}\t{miss}\n')




############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))