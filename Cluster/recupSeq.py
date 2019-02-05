#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package scriptName.py
# @author Florian CHARRIAT

"""
	The scriptName script
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
	>>> scriptName.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display scriptName.py version number and exit
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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to retrieve sequence from a list of sequence name''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-l', '--list', type = str, required=True, dest = 'list', help = 'list of sequence to process')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'path of the fasta file which contain at least all the sequence to process')
	filesreq.add_argument('-o', '--output', type=str, required=True, dest = 'output', help = 'path of output fasta file')


######### Recuperation arguments ###########
	args = parser.parse_args()
	list = os.path.abspath(args.list)
	fasta = os.path.abspath(args.fasta)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	dico = fasta2dict(fasta)
	with open(output,'w') as output_file,open(list,'r') as list_file:
		for elt in list_file :
			id = elt.split()[0].replace('>','')
			sequence = dico[id].seq
			record = SeqRecord(sequence, id=str(id), name=str(id), description= dico[id].description)
			SeqIO.write(record, output_file, "fasta")


############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in vcfSearch (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))