#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @author Florian Charriat


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens

##########Path ###################



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to retrieve the quality data of the assembly done by the ABYSS_launch script''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Quality.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-do', '--orthofinder',type = str, required=True, dest = 'pathOrthologue', help = 'Path of directory that contains all the orthofinder result (Orthogroups.GeneCount.csv,Orthogroups.txt)')
	filesreq.add_argument('-fa', '--pathFasta',type = str, required=True, dest = 'pathFasta', help = 'Path of fasta file that contains all the protein sequence of all Isolat')
	filesreq.add_argument('-o', '--outFile',type = str, required=True, dest = 'output', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	pathFasta =  os.path.abspath(args.pathFasta)
	pathOrthologue = os.path.abspath(args.pathOrthologue)
	output = os.path.abspath(args.output)



	f = open(pathOrthologue+'/Orthogroups.GeneCount.csv','r')
	lines = f.readlines()
	f.close() 
	entete = lines[0]
	listeProtein = []
	for elt in entete.split('\t') :
		if 'Total' not in elt :
			listeProtein.append(elt.replace('_protein',''))	
			
	OG_liste = lines[1:len(lines)]
	dico_OG = {}
	for line in OG_liste :
		OG = line.split('\t')[0]
		count = line.split('\t')[1:-1]
		if '0' not in count :
			dico_OG[OG] = listeProtein
		else :
			listeIndex = indexDif(count,'0')
			listeCount = []
			for i in listeIndex :
				listeCount.append(listeProtein[i])
			dico_OG[OG] = listeCount
		
	f = open(pathOrthologue+'/Orthogroups.txt','r')
	lines = f.readlines()
	f.close()
	lines = lines[1:len(lines)]
	dico_secretome = fasta2dict(pathFasta)
	dico_groupe = {}
	for line in lines :
		listeLine = line.split()
		OG = listeLine[0].replace(':','')
		if OG not in dico_OG.keys() :
			break
		Idseq = listeLine[1].replace(' ','')
		fasta_output = output +'/'+ OG + '.fasta'
		if len(dico_OG[OG]) != len(listeProtein):
			dico_groupe[OG] = []
			f = open(fasta_output,'w')
			sequence = Seq(str(dico_secretome[Idseq].seq) + '*')
			record = SeqRecord(sequence,id=str(OG),name=str(OG), description= 'length : '+ str(len(sequence)))
			SeqIO.write(record,f, "fasta")
			f.close()
	
