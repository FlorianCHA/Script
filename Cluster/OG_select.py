#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package OG_select.py
# @author Florian CHARRIAT

"""
	The OG_select script
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
	>>> OG_select.py -f file.fasta -g Orthogroups.csv -o result.fasta


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display OG_select.py version number and exit
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of the fasta which contains the sequence to process
		- \-g <path/to/orthofinder/csv/files>, --ortho <path/to/orthofinder/csv/files>
						path of the csv orthofinder result
		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file
"""
version = '0.1'

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

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to retrieved all gene with orthofinder result''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', metavar="<path/to/fasta/file>",type = str, required=True, dest = 'fasta', help = 'path of the fasta which contains the sequence to process')
	filesreq.add_argument('-db', '--database', metavar="<path/to/fasta/database/file>",type = str, required=True, dest = 'db', help = 'path of the fasta which contains all sequence used for orthofinder')
	filesreq.add_argument('-g', '--ortho', metavar="<path/to/orthofinder/csv/files>",type = str, required=True, dest = 'groups', help = 'path of the csv orthofinder result')
	filesreq.add_argument('-o', '--output', metavar="<path/to/output/file>",type=str, required=True, dest = 'output', help = 'path of the output file')

######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	ortho = os.path.abspath(args.groups)
	db = os.path.abspath(args.db)
	output = os.path.abspath(args.output)
	verifFichier(fasta)
	verifFichier(ortho)
	verifFichier(db)
########### Gestion directory ##############

############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in OG_select (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	dico = {}
	listeOG = fasta2dict(fasta).keys()
	print(form(f'Chargement des données du fichier fasta contenant ','white','bold'))
	dico_fasta = fasta2dict(db)
	print()
	with open(ortho,'r') as ortho_file, open(output,'w') as f:
		entete = ortho_file.readline()
		bar = ChargingBar('Processing :' , max=14766, suffix='%(percent)d%%')
		for line in ortho_file :
			lineSplit = line.split(':')
			OG = lineSplit[0]
			listeGene = lineSplit[1].replace('\n','').split()
			for elt in listeOG :
				if elt in listeGene :
					for elt in listeGene :
						record = SeqRecord(dico_fasta[elt].seq, id=str(elt), name=str(elt),description=dico_fasta[elt].description)
						SeqIO.write(record, f, "fasta")
					break
			bar.next()
	bar.finish()



############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))