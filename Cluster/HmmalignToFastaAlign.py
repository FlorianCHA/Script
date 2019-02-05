#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package HmmAlignToFastaAlign.py
# @author Florian CHARRIAT

"""
	The HmmAlignToFastaAlign script
	===============================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to convert Hmm alignement to fasta alignement.
	Example
	-------
	>>> HmmAlignToFastaAlign.py -f Hmm.aln -o Alignement.fasta


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display HmmAlignToFastaAlign.py version number and exit
	Input mandatory infos for running:
		- \-f <path/to/Hmm/alignement/file>, --file <path/to/Hmm/alignement/file>
						path of Hmm alignement result
		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of outout file (format fasta)

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

# Import BioPython

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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to convert Hmm alignement to fasta alignement.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'file', help = 'path of Hmm alignement result')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of outout file (format fasta)')

######### Recuperation arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	verifFichier(file)

############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("     Welcome in HmmAlignToFastaAlign (Version " + version + ")     ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	dico = {}
	with open(file,'r') as aln_file :
		for line in aln_file :
			if line[0] != '#'  and line[0:2] != '//' and line[0] != '\n':
				lineSplit = line.split()
				id = lineSplit[0]
				seq = lineSplit[1].replace('\n','')
				if id not in dico.keys():
					dico[id] = seq
				else :
					dico[id] += seq


	with open(output,'w') as ouput_file :
		for elt in dico.keys():
			record = SeqRecord(Seq(dico[elt]), id=elt, name=elt, description=f'| length : {len(dico[elt])}')
			SeqIO.write(record, ouput_file, "fasta")


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))