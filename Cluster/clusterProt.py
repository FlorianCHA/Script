#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package clusterProt.py
# @author Florian Charriat

"""
	The clusterProt script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/10/2018
	:version: 0.1

	Script description
	------------------

	This program is used to select only seq with less k % of identity.

	Example
	-------

	>>> clusterProt.py -f /homedir/user/work/file.fasta -o /homedir/user/work/result.fasta -i 90

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display clusterProt.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of the fasta file to process

		- \-i <int>, --identiity <int>
						Max identity between two sequence.

		- \-o <path/to/output/file>, --outdirPath <path/to/output/file>
						path of the output file

"""

########## Module ###############
## Python modules
import argparse, os, sys

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier,fasta2dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

if __name__ == "__main__":
	version = "0.1"

	############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__,
									 description='''This program is used to select only seq with less k % of identity.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
		'display clusterProt version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', type=str, required=True, dest='file',
						  help='Path of the fasta file to process')
	filesreq.add_argument('-i', '--identity', type=float, required=False, dest='identity',
						  help='Max identity between two sequence (default = 0.90)')
	filesreq.add_argument('-o', '--output', type=str, required=True, dest='oufile',
						  help='Path of the output file')

	######### Recuperation arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.oufile)
	identity = args.identity
	########### Gestion directory ##############
	verifFichier(file)

	########### Main ###########################
	print('\nComparaison des sequences avec ucluster\n')

	file_sort = '{}_sort.fasta'.format(file.replace('.fasta',''))
	os.system('uclust --sort {} --output {} --quiet'.format(file,file_sort))
	result_uc = '{}.uc'.format(file.replace('.fasta',''))
	os.system('uclust --input {} --uc {} --id {} --quiet'.format(file_sort, result_uc,identity))
	clusterResult =  '{}_cluster.fasta'.format(file.replace('.fasta',''))
	os.system('uclust --uc2fasta {} --input {} --output {} --types S --quiet'.format(result_uc,file_sort,clusterResult))

	print('Lecture du fichier Fasta\n')
	dico = fasta2dict(file_sort)
	liste = []
	print('Lecture du fichier uc\n')
	with open(clusterResult, 'r') as f:
		for line in f :
			if line[0] == '>':
				id_prot = line.split('|')[2].split()[0].strip()
				#evalue = float(line.split('=')[1].split('|')[0].strip())
				liste.append([id_prot])#,float(evalue)])

	#liste = sorted(liste, key=lambda sousliste: sousliste[1])

	print('Fasta in process\n')
	with open(output, 'w') as f:
		for elt in liste :
			sequence = Seq(str(dico[elt[0]].seq) + '*')
			record = SeqRecord(sequence,id=str(elt[0]),name=str(elt[0]), description= dico[elt[0]].description )
			SeqIO.write(record,f, "fasta")
	os.system('rm {} {} {}'.format(file_sort,result_uc,clusterResult))







