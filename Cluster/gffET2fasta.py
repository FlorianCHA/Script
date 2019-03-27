#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package gffET2fasta.py
# @author Florian CHARRIAT

"""
	The gffET2fasta script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to create fasta file to gff ET file (isn't a normal gff)
	Example
	-------
	>>> gffET2fasta.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display gffET2fasta.py version number and exit
	Input mandatory infos for running:
		- \-g <path/to/gff/file>, --gff <path/to/gff/file>
						path of gff ET file to used
		- \-f <path/to/fasta/file>, --fasta path/to/fasta/file>
						path of assembly file
		- \-o <path/to/output/file>, --output<path/to/output/file>
						path of output file to fasta format
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
	version =  '0.1'

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme i used to search ET in gene promoter''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', type = str, required=True, dest = 'fasta', help = 'path of assembly fasta file')
	filesreq.add_argument('-g', '--gff',type = str, required=True, dest = 'gff', help = 'path of gff ET file to used')
	filesreq.add_argument('-o', '--output',type=str, required=True, dest = 'output', help = 'path of output file to fasta format')

######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	gff = os.path.abspath(args.gff)
	output = os.path.abspath(args.output)
	verifFichier(fasta)
	verifFichier(gff)
########### Gestion directory ##############

	############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form(
		"            Welcome in vcfSearch (Version " + version + ")         ",
		type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

	########### Main #####################################

	print(form(f'Openning {fasta.split("/")[-1]} file\n'))
	dico_fasta = fasta2dict(fasta)
	print(form(f'Openning {gff.split("/")[-1]} file\n'))
	with open(gff,'r') as gff_file , open(output,'w') as output_file :
		entete = gff_file.readline()
		for line in gff_file :
			Scaffold,_,_,start,end,_,brin,_,ID = line.split()
			start = int(start)
			end = int(end)
			Scaffold = Scaffold.split('_')[-1]
			# ID = ID.split('Target=')[-1]
			ID = ID.replace('ID=','').replace(';','_').split('_Identity=')[0]
			if brin == '+' :
				seq = Seq(str(dico_fasta[Scaffold].seq)[start-1:end])
			elif brin == '-' :
				seq = Seq(str(dico_fasta[Scaffold].seq)[start - 1:end]).reverse_complement()
			else :
				print(line)
			record = SeqRecord(seq, id=str(ID), name=str(ID), description='')
			SeqIO.write(record, output_file, "fasta")

	############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))