#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Add_length_to_fasta.py
# @author Florian CHARRIAT

"""
	The Add_length_to_fasta script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1

	Script description
	------------------
	This Programme is used to add length in description to fasta file


	Example
	-------
	>>> Add_length_to_fasta.py -f file.fasta -o New_file.fasta


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Add_length_to_fasta.py version number and exit
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of the fasta file
		- \-o <path/to/fasta/file/output>, --output <path/to/fasta/file/output>
						path of the fasta file output
"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human

## BioPython
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to add length in description to fasta file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', metavar="<path/to/directory>",type = str, required=True, dest = 'fasta', help = 'path of the fasta file')
	filesreq.add_argument('-o', '--out', metavar="<path/to/directory/clumpak>",type = str, required=True, dest = 'output', help = 'path of the fasta file output')

######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	output = os.path.abspath(args.output)
	verifFichier(fasta)
########### Gestion directory ##############

############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("        Welcome in Add_length_to_fasta (Version " + version + ")   ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	fasta_dico = fasta2dict(fasta)
	with open(output,'w') as f :
		for elt in sorted(fasta_dico.keys(),key = sort_human) :
			sequence = fasta_dico[elt].seq
			record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(sequence)))
			SeqIO.write(record, f, "fasta")

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))
