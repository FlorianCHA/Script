#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package RecupInfoGFF.py
# @author Florian Charriat

########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFasta , recupId,fasta2dict,sort_human


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''TThis program is used recup information of Annotation from GFF ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display RecupInfoGFF version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the gff file')
	filesreq.add_argument('-g', '--genome',type = str, required=True, dest = 'genome', help = 'path of genome file')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'outfile', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	genome = os.path.abspath(args.genome)
	output= os.path.abspath(args.outfile)


########### Gestion directory ##############
	directory = verifDir(directory,True)

if __name__ == "__main__":

	version="0.1" 

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in RecupInfoGFF (Version " + version + ")           ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

################## Main ################################

f = open(output,'w')
f.write('%s\t%s\t%s\t%s\t%s\t%s\n'%('Strain'.center(10),'Augustus'.center(10),'Braker'.center(10),'nb Gene'.center(10),'Length  mean'.center(10),'Length (mb)'.center(10)))
for folder in os.listdir(directory) :
	if folder.endswith('_merge.gff3') :	
		gffFilePath = '%s%s'%(directory,folder)
		name = folder.replace('_merge.gff3','')
		gffFile = open(gffFilePath,'r')
		lines = gffFile.readlines()
		gffFile.close()
		nbGeneA = 0
		nbGeneB = 0
		lengthGene = 0
		lengthGenome = 0
		for line in lines :
			if line[0] != '#' :
				typeLine = line.split('\t')[2]
				typeAnnotation = line.split('\t')[1]
				posStart =  line.split('\t')[3]
				posEnd =  line.split('\t')[4]
				length = int(posEnd) - int(posStart)
			
				if typeLine == 'gene' and typeAnnotation == 'augustus_BGPI' :
						nbGeneA += 1
						lengthGene = lengthGene + length
				if typeLine == 'gene' and typeAnnotation == 'Braker_BGPI' :
						nbGeneB += 1
						lengthGene = lengthGene + length
		fastaPath = genome +'/'+name+'.fasta'
		dico_fasta = fasta2dict(fastaPath)
		for elt in dico_fasta.values():
			lengthGenome = lengthGenome +len(elt.seq)	
	
		f.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(name.center(10),str(nbGeneA).center(10),str(nbGeneB).center(10),str(nbGeneB+nbGeneA).center(10),str(round(lengthGene/(nbGeneA+nbGeneB),2)).center(10),str(round(lengthGenome/1000000,2)).center(10)))
		print('%s Done'%name)
f.close()
