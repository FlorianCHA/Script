#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package searchET.py
# @author Florian CHARRIAT

"""
	The searchET script
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
	>>> searchET.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display searchET.py version number and exit
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
	filesreq.add_argument('-f', '--fasta', type = str, required=True, dest = 'fasta', help = 'path of gene fasta file')
	filesreq.add_argument('-g', '--gff',type = str, required=True, dest = 'gff', help = 'path directory which contain all gff file ')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'directory', help = 'path of directory which contain all gff file for ET ')
	filesreq.add_argument('-o', '--output',type=str, required=True, dest = 'output', help = 'path of the output file')

######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	gff = os.path.abspath(args.gff)
	directory = os.path.abspath(args.directory)
	output = os.path.abspath(args.output)
	verifFichier(fasta)
	gff = verifDir(gff, True)
	directory = verifDir(directory,True)
########### Gestion directory ##############

	############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form(
		"            Welcome in vcfSearch (Version " + version + ")         ",
		type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

	########### Main #####################################

	dico_fasta = fasta2dict(fasta)
	dico_id = {}
	print(form(f"Ouverture du fichier {fasta}\n", 'white', 'bold'))
	for elt in sorted(dico_fasta.keys(),key = sort_human) :
		isolat = elt.split('_')[1]
		if isolat in dico_id.keys() :
			dico_id[isolat][elt] = []
		elif isolat not in dico_id.keys():
			dico_id[isolat] = {}
			dico_id[isolat][elt] = []

	with open(output,'w') as output_file :
		for elt in dico_id.keys() :
			ET_path = f'{directory}{elt}_ET.gff'
			if os.path.exists(ET_path) :
				print(form(f"Traitement de l'isolat {elt}\n",'white','bold'))

				dico_ET = {}
				print(form(f"Ouverture du fichier {ET_path}\n",'white','bold'))

				with open(ET_path,'r') as ET_file :
					for line in ET_file :
						if line[0]  !='#' :
							Scaffold, _,_, start, stop, _, brin, *_, id = line.split('\t')
							id = id .split('"')[1]
							if Scaffold not in dico_ET.keys() :
								dico_ET[Scaffold] = [(start,stop,brin,id)]
							else :
								dico_ET[Scaffold].append((start,stop,brin,id))




				gff_path = f'{gff}{elt}/{elt}_merge.gff3'
				print(form(f"Ouverture du fichier {gff_path}\n",'white','bold'))

				with open(gff_path,'r') as gff_file :
					for line in gff_file :
						if line[0] != '#' :
							Scaffold,_,type,start,stop,_,brin,*_,id = line.split('\t')
							id = id.split(';')[0].replace('ID=','')
							if type == 'mRNA' and id in dico_id[elt].keys() :
								if Scaffold in dico_ET.keys() :
									listeET = dico_ET[Scaffold]
									for ET in listeET :
										start_ET, stop_ET, brin_ET, id_ET = ET
										start = int(start)
										stop = int(stop)
										start_ET = int(start_ET)
										stop_ET = int(stop_ET)

										if brin == '+' :
											if start > start_ET > (start - 400) :
												output_file.write(f"\t - {id} have the {id_ET} element transposable dans les 400 nucléotides en amont\n")
											elif start > stop_ET> (start - 400) :
												output_file.write(f"\t - {id} have a part of the  {id_ET} element transposable dans les 400 nucléotides en amont\n")
										if brin == '-' :
											if stop < stop_ET < (stop + 400) :
												output_file.write(f"\t - {id} have the {id_ET} element transposable dans les 400 nucléotides en amont\n")
											elif stop < start_ET< (stop + 400) :
												output_file.write(f"\t - {id} have a part of the  {id_ET} element transposable dans les 400 nucléotides en amont\n")
			else :
				print(form(f"Il n'y a pas de fichier d'annotation des éléments transposable de l'isolat {elt} dans le repertoire {directory}\n", 'orange', 'bold'))











############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))