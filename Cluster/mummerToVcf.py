#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package mummerToVcf.py
# @author Florian CHARRIAT

"""
	The mummerToVcf script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to convert mummer show-snp result to vcf file
	Example
	-------
	>>> mummerToVcf.py -f show_snp_ouput -o output.vcf


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display mummerToVcf.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/show-snp/result>, --file <path/to/show-snp/result>
						path of the mummer show-snp result file
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


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to convert mummer show-snp result to vcf file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'file', help = 'path of the mummer show-snp result file')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of the output file ( format vcf)')


######### Recuperation arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	verifFichier(file)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("            Welcome in mummerToVcf (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	print(form(f'Openning show-snp result file : {file}','white','bold'))
	nb = 0
	RefFile = '/homedir/gladieux/work/magMax_project/9_SNPcalling/70-15.fasta'
	nameFile = output.split('/')[-1].split('.')[0]
	with open(file,'r') as snp_file, open(output,'w') as output_file :
		output_file.write('##fileformat=VCFv4.2\n')
		output_file.write(f'##reference=file={RefFile}\n')
		output_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		output_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{nameFile}\n')

		for line in snp_file :
			if is_number(line[0]) :
				nb += 1
				position,REF,ALT,*_,CHROM,_ = line.split()
				output_file.write(f'{CHROM}\t{position}\t.\t{REF}\t{ALT}\t.\t.\t.\tGT\t1\n')

	print(form(f'\nThe {nb} SNPs is add to {output}', 'white', 'bold'))
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))