#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package alnToSnp.py
# @author Florian CHARRIAT

"""
	The alnToSnp script
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
	>>> alnToSnp.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display alnToSnp.py version number and exit
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
from progress.bar import ChargingBar

## Python modules
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to keep only SNP from a alignement''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--alignement', type = str, required=True, dest = 'file', help = 'path of alignement file')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of output file')
	filesreq.add_argument('-t', '--type', type=str, required=True, dest = 'type', help = 'Alignement type')

######### Recuperation arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
	type =  args.type
########### Gestion directory ##############
	verifFichier(file)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in alnToSnp (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	print("\n Récupération du fichier d'alignement\n")
	alignment = AlignIO.read(file, type)
	nbSeq = len(alignment)
	length_Seq = len(alignment[1,:].seq)
	bar = ChargingBar("Traitement du fichier d'alignement", max=length_Seq, suffix='%(percent)d%%')
	nb = nbG = 0
	for i in range(0,length_Seq) :
		str_columns = alignment[:,i]
		nb += 1
		if nbSeq != str_columns.count(str_columns[0]) and 'N' not in str_columns and '-' not in str_columns:
			nbG += 1
			try :
				edited += alignment[:, i:i+1]
			except NameError:
				edited = alignment[:, i:i+1]
		bar.next()
	bar.finish()
	AlignIO.write(edited,output, 'fasta')
	print(form(f'The script kept {nbG} nucleotide of the {nb} total nucléotide of the alignement','green','bold'))
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))