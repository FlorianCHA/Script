#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package recupSeqFromGFF.py
# @author Florian Charriat

"""
	The recupSeqFromGFF script
	==========================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to recup Seq from GFF format

	Example
	-------

	>>> recupSeqFromGFF.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display recupSeqFromGFF.py version number and exit
						
	Input mandatory infos for running:
		- \-g <path/to/directory>, --gff <path/to/gff/file>
						path of gff file
		- \-s <str>, --gff<str>
						Name of strain
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
import argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from module_Flo import verifDir, createDir , form, recupId

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to recup Seq from GFF format ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display recupSeqFromGFF version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-g', '--gff',type = str, required=True, dest = 'gff', help = 'path of gff file')
	filesreq.add_argument('-s', '--strain', type=str, required=True, dest = 'strainName', help = 'Name of strain')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	gff = os.path.abspath(args.gff)
	outDir= os.path.abspath( args.outdirPath)
	strain = args.strainName

########### Gestion directory ##############
	createDir(outDir)
	outDir = verifDir(outDir)




############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in recupSeqFromGFF (Version " + version + ")        ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################

idGenome = strain
annotation = open(gff,'r')
lines = annotation.readlines()
annotation.close()
#os.system('cp augustus.gff %s.gff'%idGenome)
adn = False
protein = False 
fastaAdn = open('%s%s_transcript.fasta'%(outDir,idGenome),'w')
fastaProt = open('%s%s_protein.fasta'%(outDir,idGenome),'w')



for line in lines :	
	if 'AUGUSTUS\ttranscript\t' in line :
		ids = '%s'%line.split('\t')[8][:-1]
		numeroGene = ids.split('.t')[0].replace('g','')
		numT = str(line.split('\t')[8].split('.t')[-1]).replace('\n','')
		NewnumT = str(int(numT)-1)
		idsT = ids.replace('g'+numeroGene,"Mo_"+strain+"_"+str(numeroGene.zfill(5))+'0')
		idsT = idsT.replace('.t'+numT,'T'+NewnumT)
		idsP = ids.replace('g'+numeroGene,"Mo_"+strain+"_"+str(numeroGene.zfill(5))+'0')
		idsP = idsP.replace('.t'+numT,'P'+NewnumT)

		des = '| pos=%s:%s-%s | braker_BGPI_v1 |'%(line.split('\t')[0],line.split('\t')[3],line.split('\t')[4])
		start_gene = False
	elif '# Evidence for and against' in line :
		record = SeqRecord(Seq(seq[:-1]),id=idsT,name=idsT, description='%slength=%s'%(des,len(seq[:-1])))
		SeqIO.write(record,fastaAdn, "fasta")
		recordP = SeqRecord(Seq(prot[:-1]),id=idsP,name=idsP, description='%slength=%s'%(des,len(prot[:-1])))
		SeqIO.write(recordP,fastaProt, "fasta")
		protein = False 
	elif '# protein sequence' in line :
		prot = line.split('[')[1][:-1]
		protein = True
		adn = False
	elif '# coding sequence' in line :
		seq = line.split('[')[1][:-1]
		adn = True
	elif adn :
		seq = seq + line[2:-1]
	elif protein :
		prot = prot + line[2:-1]
	
fastaAdn.close()
fastaProt.close()

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))
