#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package OG2fasta.py.py
# @author Florian CHARRIAT

"""
	The OG2fasta.py script
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
	>>> OG2fasta.py.py -f multiSample.vcf -o output.vcf -d show-snp.result


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display OG2fasta.py.py version number and exit

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
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,is_number,isIn

## Python modules
import argparse
from time import gmtime, strftime
import pandas as pd

## BIO Python modules
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to correct vcf file with multi sample create with show-snp''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'file', help = 'path of thec fasta file wich contains all sequence to group by OG')
	filesreq.add_argument('-g', '--group',type = str, required=True, dest = 'group', help = 'path of the output file of Orthofinder (Orthogroups.txt)')
	filesreq.add_argument('-db', '--database', type=str, required=True, dest='db',
						  help='path of the database file which contain all seq at cds format')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of the output directory')


######### Recuperation arguments ###########
	args = parser.parse_args()
	group = os.path.abspath(args.group)
	db = os.path.abspath(args.db)
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	verifFichier(file)
	verifFichier(group)
	directory = verifDir(output)
	createDir([output,f'{output}/fasta'])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("           Welcome in OG2fasta.py (Version " + version + ")        ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	print('Lecture fichier Orthogroups')
	dico_fasta = fasta2dict(file)
	listeSeq = dico_fasta.keys()
	dico_OG = {}
	with open(group, 'r') as group_file :
		for line in group_file :
			OG,gene = line.split()[0].replace(':',''), line.split()[1:len(line.split())]
			if isIn(gene,listeSeq) :
				dico_OG[OG] = gene

	print('Récupération de la DB')
	db_dico = fasta2dict(db)
	print('Récupération séquence')
	liste = []
	for OG in dico_OG.keys() :
		with open(f'{output}/fasta/{OG}.fasta','w') as output_file :
			for seq in dico_OG[OG] :
				name = seq
				if 'Mo_' in seq :
					name = seq.replace('Mo_','')
				name = name.split('_')[0]
				liste.append(name)
				record = SeqRecord(seq= db_dico[seq].seq, id=name, name=name, description='')
				SeqIO.write(record, output_file, "fasta")

	liste = set(liste)
	with open(f'{output}/group.txt', 'w') as output_file:
		for elt in sorted(liste,key = sort_human) :
			output_file.write(f'{elt}\n')
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))