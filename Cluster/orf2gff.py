#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package orf2gff.py
# @author Florian Charriat

"""
	The orf2gff script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 15/10/2018
	:version: 0.1

	Script description
	------------------

	This program is used to create a gff file from a getORF result and the fasta file.

	Example
	-------

	>>> orf2gff.py -f /homedir/user/work/fasta.file -p /homedir/user/work/getORF.result -o /homedir/user/work/result.gff

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display orf2gff.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta that contains all the sequence to use.
						
		- \-p <path/to/HMM/profil/file>, --profil <path/to/HMM/profil/file>
						path of the getORF result which be used

		- \-o <path/to/output/directory>, --output <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form,verifFichier,fasta2dict, openfile



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to create a gff file from a getORF result and the fasta file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display orf2gff version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of fasta that contains all the sequence to use.')
	filesreq.add_argument('-orf', '--getORF',type = str, required=True, dest = 'orf', help = 'Path of the getORF result which be used')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	orf = os.path.abspath(args.orf)
	output = os.path.abspath( args.output)

########### Gestion directory ##############
	verifFichier(fasta)
	verifFichier(orf)
	os.system('mkdir %s'%output)

################# Main #####################

	listeORF = openfile(orf)
	fasta = openfile(fasta)
	ids = ''
	startfile = True
	dico = {}
	dico_brin = {}
	listeId = []
	for sequence in fasta:
		if sequence[0] == '>' :
			listeId.append(sequence.split()[0])
			dico_brin[sequence.split()[0].split(':N')[0]] = sequence.split('brin :')[1].split(',')[0]

	def notIn(elt1,elt2):
		for elt in elt2 :
			if elt[0] <= elt1[0] <= elt[1] or elt[0] <= elt1[1] <= elt[1] or elt1[0] <= elt[1] <= elt1[1] or elt1[0] <= elt[1] <= elt1[1]:
				return False
		
		return True

	def bestcouverture(liste):
		listeProvisoire = []
		for elt in liste :	
			elt[0] = int(elt[0])
			elt[1] = int(elt[1])
			length = elt[1]-elt[0]
			listeProvisoire.append([elt[0],elt[1],length])
		listeProvisoire = sorted(listeProvisoire ,key=lambda colonnes: colonnes[2], reverse = True)
		listeFinal = [listeProvisoire[0]]
		for elt in listeProvisoire :
			if notIn(elt,listeFinal) :
				listeFinal.append(elt)
		return listeFinal

	i = 0
	for orf in listeORF :
		if orf[0] == '>' and '(REVERSE SENSE)' not in orf:
			orfSplit = orf.split('[')[1].split(']')[0].split('-')
			start = orfSplit[0]
			end = orfSplit[1]
			if ids != orf.split('_')[0].replace('>','') and startfile == False:
				i += 1
				ids = orf.split('_')[0].replace('>','')
				if listeId[i].split(':')[-1] == ids :
					idFasta = listeId[i].split(':N')[0]
					dico[idFasta] = [[start,end]]
				else :
					print(listeId[i].split(':')[-1],ids)
			elif startfile == True:
				ids = orf.split('_')[0].replace('>','')
				if listeId[i].split(':')[-1] == ids :
					idFasta = listeId[i].replace('>','').split(':N')[0]
					dico[idFasta] = [[start,end]]
				else :
					print(listeId[i].split(':')[-1],ids)
			else :
				ids = orf.split('_')[0].replace('>','')
				dico[idFasta].append([start,end])
			startfile = False



for ids in dico.keys() :
	cds = bestcouverture(dico[ids])
	isolat = ids.split('_')[0].replace('>','')
	f = open(output+'/'+isolat+'.gff','a')
	if '>' not in ids :
		brin = dico_brin['>'+ids].strip()
	else :
		brin = dico_brin[ids].strip()
	start = int(ids.split(':')[1])
	stop = int(ids.split(':')[2])
	Scaffold = 'Scaffold_%s'%ids.split('_')[2].split(':')[0]
	f.write('%s\tPutative\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n'%(Scaffold,start,stop,brin,ids,ids))
	f.write('%s\tPutative\tmRNA\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start,stop,brin,ids,ids))
	if brin== '+' :
		f.write('%s\tPutative\tstart_codon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start,start+2,'.',ids,ids))
		for elt in cds :
			debut = elt[0] - 1
			fin = elt[1] - 1
			f.write('%s\tPutative\tCDS\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start+debut,start+fin,brin,ids,ids))
			f.write('%s\tPutative\texon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start+debut,start+fin,brin,ids,ids))
		f.write('%s\tPutative\tstop_codon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,stop-2,stop,'.',ids,ids))

	if brin == '-' :
		f.write('%s\tPutative\tstop_codon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start,start+2,'.',ids,ids))
		for elt in cds :
			debut = elt[0] - 1
			fin = elt[1] - 1
			f.write('%s\tPutative\tCDS\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start+debut,start+fin,brin,ids,ids))
			f.write('%s\tPutative\texon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,start+debut,start+fin,brin,ids,ids))
		f.write('%s\tPutative\tstart_codon\t%s\t%s\t.\t%s\t.\tID=%sT0;Parent=%s\n'%(Scaffold,stop-2,stop,'.',ids,ids))










		





