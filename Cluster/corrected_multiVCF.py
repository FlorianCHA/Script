#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package corrected_multiVCF.py
# @author Florian CHARRIAT

"""
	The corrected_multiVCF script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to correct vcf file with multi sample create with show-snp
	Example
	-------
	>>> corrected_multiVCF.py -f multiSample.vcf -o output.vcf -d show-snp.result


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display corrected_multiVCF.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/show-snp/result>, --file <path/to/show-snp/result>
						path of the vcf file with multi sample
		- \-d <path/to/show-snp/directory>, --directory <path/to/show-snp/directory>
						path of the directory which contains all show-snp file without -I and -C option.
		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file ( format vcf)

"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,is_number

## Python modules
import argparse
import pandas as pd


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to correct vcf file with multi sample create with show-snp''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'file', help = 'path of the vcf file with multi sample')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'directory', help = 'path of the directory which contains all show-snp file without -I and -C option.')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of the output file ( format vcf)')


######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.directory)
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	verifFichier(file)
	verifDir(directory,True)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in corrected_multiVCF (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	print(form(f'Openning vcf file with multi sample : {file}','white','bold'))
	test = pd.read_csv(file, header=None)
	print(test[1:])



	print(form(f'\nThe {nb} SNPs is add to {output}', 'white', 'bold'))
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))