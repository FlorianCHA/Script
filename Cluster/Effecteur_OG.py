#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Effecteur_OG.py
# @author Florian CHARRIAT

"""
	The Effecteur_OG script
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
	>>> Effecteur_OG.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Effecteur_OG.py version number and exit
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
import sys, os,re
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,isIn

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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to analysed effecteur group orthologue''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of effecteur fasta file')
	filesreq.add_argument('-o', '--orthogroup', type = str, required=True, dest = 'group', help = 'Path of the result of orthofinder (format txt)')
	filesreq.add_argument('-db', '--db',type=str, required=True, dest = 'db', help = 'Path to fasta file which be used for orthofinder analysis')
	filesreq.add_argument('--output',type=str, required=True, dest = 'output', help = 'Path to output directory')


######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	gr = os.path.abspath(args.group)
	db = os.path.abspath(args.db)
	output = os.path.abspath(args.output)

########### Gestion directory ##############
	verifFichier(fasta)
	verifFichier(gr)
	verifFichier(db)
	output = verifDir(output)
	createDir([output])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in Effecteur_OG (Version " + version + ")       ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	dico_effecteur = fasta2dict(fasta)
	entete = []
	for elt in dico_effecteur :
		isolat = elt.split('_')[1]
		if isolat not in entete :
			entete.append(isolat)

	effecteur = dico_effecteur.keys()
	dico = {}
	print('Création du fichier Effecteur_count')
	with open(gr,'r') as orthofinder_result, open(f'{output}Effecteur.count.csv','w') as count_file :
		count_file.write('\t'.join(entete)+"\n")
		for line in orthofinder_result :
			OG = line.split(':')[0]
			genes = line.split()[1:len(line.split())]
			line_output = f'{OG}'
			if isIn(effecteur, genes):
				dico[OG] = genes
				for elt in sorted(entete,key = sort_human) :
					pattern = f'[^" "]*{elt}[^" "]*'
					nb_isolat = len(re.findall(pattern, line))
					line_output = line_output + f'\t{nb_isolat}'
				count_file.write(line_output +'\n')
	print('\nCréation du fichier des distribution de longueur pour chaque groupe Orthologue')
	fasta_bd = fasta2dict(db)
	with open(f'{output}Effecteur_length.txt', 'w') as len_file:
		len_file.write(f'Orthogroup\tlength\n')
		for elt in dico.keys():
			for gene in dico[elt] :
				len_file.write(f'{elt}\t{len(fasta_bd[gene])}\n')


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))