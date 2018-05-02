#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package eliminateREmotif.py
# @author Florian Charriat

"""
	The eliminateREmotif script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/04/2018
	:version: 0.1

	Script description
	------------------

	This program is used to eliminate the result of PS-scan with PS00014 pattern (RE target Motif). This program is used by secretome_Pipeline.
	
	Example
	-------

	>>> eliminateREmotif.py  -f /homedir/user/work/data/TMHMM_result.txt -o /homedir/user/work/result/


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display selection_TMHMM.py version number and exit
						
	Input mandatory infos for running:
								
		- \-p <path/to/ps_scan/output/file>, --ps_scan <path/to/ps_scan/output/file>
						path of the ps_scan output file
		- \-f <path/to/fasta/file>, --file <path/to/fasta/file>
						Path of the fasta file which ps_scan has been proceed
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path and name of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''
	This program is used to eliminate the result of PS-scan with PS00014 pattern (RE target Motif). This program is used by secretome_Pipeline.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display eliminateREmotif version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-p', '--ps_scan',type = str, default = 'None', dest = 'ps_scan', help = 'Path of the ps_scan output file')
	filesreq.add_argument('-f', '--fasta',type = str,  required=True, dest = 'fasta', help = 'Path of the fasta file which ps_scan has been proceed')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdir', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	outDir= os.path.abspath(args.outdir)
	PScan = os.path.abspath(args.ps_scan)
	fasta = os.path.abspath(args.fasta)
	Id = recupId(fasta.split('/')[-1])
	outDir = verifDir(outDir,check = False)
########### Gestion directory ##############

	dico = fasta2dict(fasta)
	f = open(PScan,'r')
	lines = f.readlines()
	f.close()
	for line in lines :
		line = line.split('\t')
		del dico[line[0]]
		
	f = open('%s%s_secreted_3.fasta'%(outDir,Id),'w')
	for idSeq in sorted(dico.keys(), key=sort_human):
			seqObj = dico[idSeq].seq
			record = SeqRecord(seqObj,id=idSeq,name=idSeq, description= dico[idSeq].description)
			SeqIO.write(record,f, "fasta")

	
	f.close()

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in selection_TMHMM   (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')
	
	
	
