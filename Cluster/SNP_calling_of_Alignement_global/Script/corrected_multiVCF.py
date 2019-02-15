#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package corrected_multiVCF.py
# @author Florian CHARRIAT

"""
	The corrected_multiVCF script
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
	>>> corrected_multiVCF.py -f multiSample.vcf -o output.vcf -d show-snp.result


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display corrected_multiVCF.py version number and exit

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
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,is_number

## Python modules
import argparse
from time import gmtime, strftime
import pandas as pd


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
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'file', help = 'path of the vcf file with multi sample')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'directory', help = 'path of the directory which contains all show-snp file without -I and -C option.')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of the output file ( format vcf)')


######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.directory)
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	verifFichier(file)
	directory = verifDir(directory,True)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("       Welcome in corrected_multiVCF (Version " + version + ")     ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	print(form(f'Openning vcf file with multi sample : {file}','white','bold'))

	with open(file,'r') as multi_file, open(output,'w') as output_file:
		for line in multi_file :
			if line[0:2] == '##' :
				output_file.write(line)
			if line[0:2] == '#C' :
				colnames = line.replace('#','').replace('\n','').split('\t')
				output_file.write('#')
				break
	multivcf = pd.read_csv(file, header=15,sep = '\t',names = colnames)
	for elt in range(9,multivcf.shape[1]) :
		name = multivcf.columns[elt]
		dico_isolat = {}
		print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

		print(form(f"Traitement de l'échantillon {name}", 'Green', 'bold'))
		print()
		print(form(f"\t- Ouverture du fichier {name}_all.snp", 'white', 'bold'))
		with open(f'{directory}{name}_all.snp','r') as all_snp :
			for line in all_snp :
				p1,ref,alt,P2,*_,K,_ = line.split('\t')
				if K not in dico_isolat.keys() :
					dico_isolat[K] = {}
				dico_isolat[K][int(p1)] = alt
		print(form(f"\t- Correction des données de l'échantillon {name}", 'white', 'bold'))
		for i in range(0,multivcf.shape[0]):
			if multivcf.iloc[i,elt] == '.' :
				K = multivcf.iloc[i,0]
				position = int(multivcf.iloc[i,1])
				if position in dico_isolat[K].keys() :
					if dico_isolat[K][position] == '.' :
						multivcf.iat[i, elt] = '.'
					if dico_isolat[K][position] != '.' :
						multivcf.iat[i, elt] = 'N'
				else :
					multivcf.iat[i, elt] = '0'

	multivcf.to_csv(output,sep='\t',mode = 'a', index = False)
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))