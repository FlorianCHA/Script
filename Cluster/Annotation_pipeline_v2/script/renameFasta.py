#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package renameFasta.py
# @author Florian Charriat

"""
	The renameFasta script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 04/07/2018
	:version: 0.1

	Script description
	------------------

	This program is used to rename output of gff2fasta.pl for the pipeline Annotation_pipeline.snake
	
	Example
	-------

	>>> renameFasta.py  -d /work/outputDirectory/ -s stain


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display renameFasta.py version number and exit
						
	Input mandatory infos for running:
								
		- \-d <path/to/gff2fasta.pl/output/directory>, --file <path/to/augustus/output/file>
						path of the gff2fasta output directory
		- \-s <strain/name>, --strain <strain/name>
						Name of strain



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
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to rename output of gff2fasta.pl for the pipeline Annotation_pipeline.snake''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display renameFasta version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str ,  required=True,dest = 'dir', help = 'Path of the gff2fasta.pl output directory')
	filesreq.add_argument('-s', '--strain',type = str ,  required=True,dest = 'strain', help = 'Name of the strain')


	
######### Recuperation arguments ###########

	args = parser.parse_args()
	directory = os.path.abspath(args.dir)
	folder= args.strain
############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in renameFasta   (Version " + version + ")         ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############## Main ##################################

	dico_Gff = {}
	cdsFile = directory +"/"+ folder + '_cds.fna'
	protFile = directory +"/"+ folder + '_protein.faa'
	gffFile = directory +"/"+ folder + '_merge.gff3'
	f = open(gffFile,'r')
	lines = f.readlines()
	f.close()
	print('Creation dico GFF')
	for line in lines :
		if 'mRNA' in line :
			lineSplit = line.split('\t')
			ids = lineSplit[8].split('ID=')[-1].split(';Parent=')[0]
			position =  'pos=%s_%s:%s'%(lineSplit[0],lineSplit[3],lineSplit[4])
			tools = lineSplit[1]
			dico_Gff[ids] = (position,tools)
	print('Creation dico des fasta')
	dico_cds = fasta2dict(cdsFile)
	dico_prot = fasta2dict(protFile)
	f = open(cdsFile.replace('.fna','.fasta'),'w')
	for idSeq in sorted(dico_cds.keys(), key=sort_human):
		position = dico_Gff[idSeq][0]
		tools = dico_Gff[idSeq][1]
		length = len(str(dico_cds[idSeq].seq))
		seqObj = dico_cds[idSeq].seq
		record = SeqRecord(seqObj,id=idSeq,name=idSeq, description='| %s | %s | %s ' %(position,tools,length))
		SeqIO.write(record,f, "fasta")
	f.close()
	f = open(protFile.replace('.faa','.fasta'),'w')
	for idSeq in sorted(dico_prot.keys(), key=sort_human):
		position = dico_Gff[idSeq.replace('0P','0T')][0]
		tools = dico_Gff[idSeq.replace('0P','0T')][1]
		length = len(str(dico_prot[idSeq].seq))
		seqObj = Seq(str(dico_prot[idSeq].seq).replace('*',''))
		record = SeqRecord(seqObj,id=idSeq,name=idSeq, description='| %s | %s | length=%s ' %(position,tools,length))
		SeqIO.write(record,f, "fasta")
	f.close()
