#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package clusterProt.py
# @author Florian Charriat

"""
	The clusterProt script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/10/2018
	:version: 0.1

	Script description
	------------------

	This program is used to select only seq with less k % of identity.

	Example
	-------

	>>> clusterProt.py -f /homedir/user/work/file.fasta -o /homedir/user/work/result.fasta -i 90

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display clusterProt.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of the fasta file to process

		- \-i <int>, --identiity <int>
						Max identity between two sequence (default = 0.90).

		- \-o <path/to/output/file>, --outdirPath <path/to/output/file>
						path of the output file

"""

########## Module ###############
## Python modules
import argparse, os, sys

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier,fasta2dict,sort_human
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os


def filtreSequence(fasta,file_out,mean,intervalle):
	"""
	Remove only the redundante sequence
	"""
	dico = fasta
	listeSequence = []
	nb = 0
	nbT = 0
	nbE = 0
	with open(file_out,'w') as f :
		for elt in sorted(dico.keys(),key = sort_human) :
			seq = str(dico[elt].seq)
			nbT += 1
			if seq not in listeSequence :
				nbE +=1
				if intervalle > 0 and (mean - intervalle) < len(seq) < (mean + intervalle) or intervalle < 0 :
					nb += 1
					listeSequence.append(seq)
					sequence = Seq(seq)
					record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(seq)))
					SeqIO.write(record, f, "fasta")

	return (nb/nbT)*100, (1 - (nbE/nbT))*100

if __name__ == "__main__":
	version = "0.1"

	############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__,
									 description='''This program is used to select only seq with less k % of identity.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
		'display clusterProt version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', type=str, required=True, dest='file',
						  help='Path of the fasta file to process')
	filesreq.add_argument('-o', '--output', type=str, required=True, dest='oufile',
						  help='Path of the output file')

	files = parser.add_argument_group('Input infos for running with default values')
	filesreq.add_argument('-i', '--identity', type=float, required=False, dest='identity',
						  help='Max identity between two sequence (default = 0.90)')
	filesreq.add_argument('--noredundant',action='store_true', dest='noredundant',
						  help='Remove only redundant sequence')
	files.add_argument('-q', '--quiet',action='store_true',
					   dest='quiet', help='No stdout if you use this option')
	files.add_argument('--intervalle', type=int,required=False, default = -1,
					   dest='inter', help='length of sequence')
	######### Recuperation arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.oufile)
	identity = args.identity
	noredundant = args.noredundant
	inter = args.inter

	quiet = args.quiet
	########### Gestion directory ##############
	verifFichier(file)

	########### Main ###########################

	if noredundant and inter > 0 :
		if quiet == False:
			print('\nEliminate only redundant sequence\n')
		fasta = fasta2dict(file)
		liste_len = []
		for elt in fasta :
			liste_len.append(len(fasta[elt].seq))
		mean_len = round((sum(liste_len) / len(liste_len)),0)
		print(f'Length mean : {mean_len}')
		nb,nbe = filtreSequence(fasta,output,mean_len,inter)
		nb = round(nb,2)
		nb = round(nb, 2)
		print(f'This script kept {nb} % of sequence ( {nbe} % redundant)')
		exit()



	elif noredundant :
		listeSequence = []
		if quiet == False:
			print('\nEliminate only redundant sequence\n')
		fasta = fasta2dict(file)
		with open(output, 'w') as f:
			for elt in sorted(fasta.keys(), key=sort_human):
				seq = str(fasta[elt].seq)
				if seq not in listeSequence:
					listeSequence.append(seq)
					sequence = Seq(seq)
					record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(seq)))
					SeqIO.write(record, f, "fasta")
		exit()





	if quiet == False :
		print('\nComparaison des sequences avec ucluster\n')
	file_sort = '{}_sort.fasta'.format(file.replace('.fasta',''))
	os.system('uclust --sort {} --output {} --quiet'.format(file,file_sort))
	result_uc = '{}.uc'.format(file.replace('.fasta',''))
	os.system('uclust --input {} --uc {} --id {} --quiet'.format(file_sort, result_uc,identity))
	clusterResult =  '{}_cluster.fasta'.format(file.replace('.fasta',''))
	os.system('uclust --uc2fasta {} --input {} --output {} --types S --quiet'.format(result_uc,file_sort,clusterResult))
	if quiet == False:
		print('Lecture du fichier Fasta\n')
	dico = fasta2dict(file_sort)
	liste = []
	if quiet == False:
		print('Lecture du fichier uc\n')
	with open(clusterResult, 'r') as f:
		for line in f :
			if line[0] == '>':
				id_prot = line.split('|')[2].split()[0].strip()
				#evalue = float(line.split('=')[1].split('|')[0].strip())
				liste.append([id_prot])#,float(evalue)])

	#liste = sorted(liste, key=lambda sousliste: sousliste[1])
	if quiet == False:
		print('Fasta in process\n')
	with open(output, 'w') as f:
		for elt in liste :
			sequence = Seq(str(dico[elt[0]].seq) + '*')
			record = SeqRecord(sequence,id=str(elt[0]),name=str(elt[0]), description= dico[elt[0]].description )
			SeqIO.write(record,f, "fasta")
	os.system('rm {} {} {}'.format(file_sort,result_uc,clusterResult))







