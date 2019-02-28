#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package makeEffectome.py
# @author Florian Charriat

"""
	The makeEffectome script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 15/10/2018
	:version: 0.1

	Script description
	------------------

	This program is used to create a fasta file which contain all nucléotidique sequence of effecteur MAX with 100bp before and after the sequence.

	Example
	-------

	>>> makeEffectome.py -f /homedir/user/work/fasta.file -p /homedir/user/work/getORF.result -o /homedir/user/work/result.gff

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display makeEffectome.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta that contains all the genome assembly of Isolat
						
		- \-p <path/to/fasta/file>, --prot <path/to/fasta/file>
						path of fasta that contains all the protein of Isolat.

		- \-hmm <path/to/hmmsearch/file>, --hmmsearch <path/to/fasta/file>
						path of hmmsearch result (stdout of effecteur MAX search)

		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file

"""


########## Module ###############
## Python modules
import argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from progress.bar import ChargingBar
#Import module_Flo
from module_Flo import verifDir, createDir , form,verifFichier,fasta2dict, openfile,is_number



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to create a fasta file which contain all nucléotidique sequence of effecteur MAX with 100bp before and after the sequence.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display makeEffectome version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of directory that contains all the genome assembly of Isolat')
	filesreq.add_argument('-g', '--gff',type = str, required=True, dest = 'gff', help = 'Path of directory that contains all the gff file of Isolat.')
	filesreq.add_argument('-list', '--list',type = str, required=True, dest = 'hmm', help = 'Path of the file witch contains all sequence ton analyse')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	genome = os.path.abspath(args.fasta)
	gff = os.path.abspath(args.gff)
	hmm = os.path.abspath(args.hmm)
	output = os.path.abspath( args.output)

########### Gestion directory ##############
	verifDir(genome,True)
	verifDir(gff,True)
	verifFichier(hmm)

################## Main ######################

	with open(hmm,'r') as list_file :
		list = list_file.read()
		list = list.replace('\n',' ').replace('\t',' ')
		list = list.split()

	for gene in list :
		if 'Mo_' not in gene :
			print(form(f"The gene {gene} haven't a good name, please rename you gene with Organisme_isolat_Nam",'red','bold'))
			exit()
		print(f'Traitement du gène : {gene}')
		isolat = gene.split('_')[1]
		fasta = fasta2dict(f'{genome}/{isolat}.fasta')
		with open(f'{gff}/{isolat}_merge.gff3','r') as gff_file, open(output,'a') as output_file:
			for line in gff_file :
				if gene in line and 'mRNA' in line :
					Scaffold,_,_,start,end,_,brin,*_ = line.split()
					start = int(start)
					end = int(end)
					if brin == '+' :
						seq = Seq(str(fasta[Scaffold].seq)[start-200:start])
					if brin == '-':
						seq = Seq(str(fasta[Scaffold].seq)[end:end +200])
						seq  =seq.reverse_complement()

					record = SeqRecord(seq, id=gene, name=gene,description ='')
					SeqIO.write(record, output_file, "fasta")

