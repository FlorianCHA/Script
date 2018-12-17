#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package slice-fasta-alignement.py
# @author Florian CHARRIAT

"""
	The slice-fasta-alignement
	===========================
	:author: Sebastien Ravel
	:contact: florian.charriat@inra.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------
	This Programme is used to slice a alignement between position 1 and position 2

	Example
	-------

	>>> slice-fasta-alignement.py -f /homedir/input.fasta -pos1 5 -pos2 100 -o /homedir/result.fasta

	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display slice-fasta-alignement.py version number and exit
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta file to process
		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output fasta file

	Input infos for running with default values:
		- \-pos1 <int>, --position1 <int>
						The first position from which the alignement is kept (default = start of alignement)
		- \-pos2 <int>, --position2 <int>
						The last position from which the alignement is kept (default = end of alignement)


"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human

## Python modules
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='27/01/2017'


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to slice a alignement between position 1 and position 2.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', metavar="<path/to/directory>",type = str, required=True, dest = 'fasta', help = 'Path of fasta file to process')
	filesreq.add_argument('-o', '--output', metavar="<path/to/output/file>",type = str, required=True, dest = 'output', help = 'Path of the output fasta file')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-pos1', '--position1', type=int, required=False,default = 1, dest = 'pos1', help = 'The first position from which the alignement is kept (default = start of alignement (position 1))')
	files.add_argument('-pos2', '--position2', type=int, required=False, default = 0,dest = 'pos2', help = 'The last position from which the alignement is kept (default = end of alignement)')
	files.add_argument('-k', '--keeporder', action='store_true', dest='k',
					   help='Keep the origin order of sequenc, otherwise a classic order is apply')

######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	output = os.path.abspath(args.output)
	pos1 = args.pos1
	pos2 = args.pos2
	keep = args.k


	########### Gestion directory ##############
	verifFichier(fasta)

############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in vcfSearch (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	aln = fasta2dict(fasta)
	if keep :
		print(form('Keep option is used. A variable is create for keep the order of sequence.','white','bold'))
		listeSequence = []
		with open(fasta,'r') as f :
			for line in f :
				if line[0] == '>' :
					listeSequence.append(line.split()[0].replace('>',''))
	else :
		listeSequence = sorted(aln.keys(), key = sort_human)
	with open(output,'w') as f :
		for elt in listeSequence :
			if pos2 == 0 :
				pos2 = len(aln[elt].seq)
			sequence = str(aln[elt].seq)[pos1-1:pos2]
			sequence = Seq(sequence)
			record = SeqRecord(sequence, id=str(elt), name=str(elt), description= aln[elt].description)
			SeqIO.write(record, f, "fasta")

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))