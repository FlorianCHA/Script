#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package mergeBraker_augustus.py
# @author Florian Charriat

"""
	The mergeBraker_augustus script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 11/06/2018
	:version: 0.1

	Script description
	------------------

	This program is used to merged Augustus ouput with braker output. This program is used in Braker_pipeline. Augustus output must be at gff3 format.
	
	Example
	-------

	>>> mergeBraker_augustus.py  -a /work/augustus.gff3 -b /work/braker.gff3 -outDir /work/result/


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display mergeBraker_augustus.py version number and exit
						
	Input mandatory infos for running:
								
		- \-a <path/to/augustus/output/file>, --file <path/to/augustus/output/file>
						path of the augustus output file
		- \-b <path/to/braker/output/file>, --file <path/to/braker/output/file>
						path of the braker output file
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
		This program is used to merged Augustus ouput with braker output. This program is used in Braker_pipeline. Augustus output must be at gff3 format.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display mergeBraker_augustus version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-a', '--augustus',type = str, default = 'None', dest = 'augustus', help = 'Path of the augustus output file')
	filesreq.add_argument('-b', '--braker',type = str,  required=True, dest = 'braker', help = 'Path of the braker output file')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdir', help = 'Path of the output directory')

	
######### Recuperation arguments ###########

	args = parser.parse_args()
	brakerFile = os.path.abspath(args.braker)
	augustusFile = os.path.abspath(args.augustus)
	outDir= os.path.abspath(args.outdir)	
########### Gestion directory ##############

	name_directory = [outDir]
	createDir(name_directory)
	outDir = verifDir(outDir)
	name = augustusFile.split('/')[-1].replace('.gff3','')
	


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in mergeBraker_augustus   (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

####################### main #################"
	nouveauGene = []
	listeGene = []
	liste = []
	liste1 = []
	nbAjout = 0
	if nbAjout != 0 :
		print('Nombre de gene identique : '+str(nbAjout))
	nbAjout = 0
	braker = open(brakerFile,'r')
	brakerlines = braker.readlines()
	braker.close()
	augustus = open(augustusFile,'r')
	augustuslines = augustus.readlines()
	augustus.close()
	print(name +' in process')
	for lineA in augustuslines :
		if lineA[0] != '#' :
			listeGene.append(lineA)
		elif 'end gene' in lineA :
			liste.append(listeGene)
			liste1.append(listeGene)
			listeGene = []
	print('Liste Augustuse done')
	print('Nombre de gene dans le fichier augustus : '+ str(len(liste)))
	for elt in liste :
		gene = elt[0]
		scaffoldA = gene.split('\t')[0]
		for lineB in brakerlines :
			scaffoldB = lineB.split('\t')[0]
			if scaffoldA == scaffoldB :
				scaffold = 'done'
				typeB = lineB.split('\t')[2]
				if typeB == 'gene' :
					startA = int(gene.split('\t')[3])
					endA = int(gene.split('\t')[4])
					startB = int(lineB.split('\t')[3])
					endB = int(lineB.split('\t')[4])
					id = lineB.split('\t')[-1].replace('\n','')
					if (startB <= endA <= endB) or (startB <= startA <= endB) or (startA <= endB <= endA) or (startA <= startB <= endA): 
						#print('>%s\t%s\t%s\t%s\t%s:%s\t%s:%s'%(id,scaffoldA,scaffoldB,typeB,startA,endA,startB,endB))
						nbAjout +=1
						liste1.remove(elt)
						break
	numberGeneBraker = brakerlines[-1].split('=g')[-1].split('.t')[0]
	print('Comparaison done')
	os.system('./renameGFF.py -g %s -s %s -o %s%s_braker.gff3'%(brakerFile,name,outDir,name))
	print('rename Braker file done')								
	print('Nombre de genes a ajouter : '+ str(len(liste1)))
	f = open(outDir+name+'_prov.gff3','w')
	for elt in liste1:
		for sousElt in elt :
			f.write(sousElt)
	f.close()
	print('Write new gff3 for augustus : Done')
	os.system('./renameGFF.py -g %s -s %s -o %s%s_augustus.gff3 -n %s -t augustus'%(outDir+name+'_prov.gff3',name,outDir,name,numberGeneBraker))
	os.system('rm %s'%(outDir+name+'_prov.gff3'))
	os.system('cat %s%s_braker.gff3 %s%s_augustus.gff3 > %s%s_merge.gff3'%(outDir,name,outDir,name,outDir,name))
			
				
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




			
							
				
				
				
				
			
				
				
				
				
				
							
					
									
									
