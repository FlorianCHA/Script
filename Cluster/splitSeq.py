#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package splitSeq.py
# @author Florian CHARRIAT

"""
	The splitSeq script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/07/2018
	:version: 0.1

	Script description
	------------------

	This program is used to corrected some Orthofinder probleme (2 Genes in a same sequence, so the orthogroups contains 2 sequences)
	
	Example
	-------

	>>> splitSeq.py -f /homedir/user/work/file.fasta -o /homedir/user/work/result/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to//fasta/file/>, --fasta <path/to//fasta/file/>
						Path of the fasta which contains one orthogroups to corrected
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from io import StringIO
from Bio import AlignIO
## Python modules
import argparse


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	version = "0.1"

############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to corrected some Orthofinder probleme (2 Genes in a same sequence, so the orthogroups contains 2 sequences''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= 'display script_pangenome version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta', type=str, required=True, dest='fasta', help='Path of the fasta which contains one orthogroups to corrected')
	filesreq.add_argument('-p', '--prefix', type=str, required=True, dest='prefix', help='Prefix name')

	######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	verifFichier(fasta)
	prefix =  os.path.abspath(args.prefix)
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print("\t" + form("|", 'yellow', 'bold') + form("            Welcome in splitSeq (Version " + version + ")          ",type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	try :
		clsutal_aln = ClustalwCommandline("clustalo",infile=fasta, outfile = fasta.replace('.fasta','.aln'))
		stdout, stderr = clsutal_aln()
	except :
		print(form('This script use the Clustalo tools, please charge this tools (for the Cirad cluster use : module load bioinfo/clustalo/1.2.0)', 'red','bold'))

	align = AlignIO.read(fasta.replace('.fasta','.aln'), "clustal")
	# align = AlignIO.read(StringIO(stdout), "clustal")
	dico_aln_rigth = {} # Dico contenant les alignements pour les quels les gaps situé avant la séquence seront éliminer
	dico_aln_left = {} # Dico contenant les alignements pour les quels les gaps situé aprés la séquence seront éliminer
	dico = {}
	for elt in align :
		dico_aln_rigth[elt.id] = str(elt.seq).lstrip('-')
		dico_aln_left[elt.id] = str(elt.seq).rstrip('-')
		dico[elt.id] = str(elt.seq)

	min_rigth =  len(align[1].seq) # Récupéré la longueurs de séquences la plus grande pour chaque min
	min_left = len(align[1].seq)
	for elt in dico_aln_rigth.keys() :
		if len(dico_aln_rigth[elt]) < min_rigth :
			min_rigth = len(dico_aln_rigth[elt])
			id_min_rigth = elt
		if len(dico_aln_left[elt]) < min_left :
			min_left = len(dico_aln_left[elt])
			id_min_left = elt

	seq_min_rigth = dico_aln_rigth [id_min_left]
	len_rigth = len(seq_min_rigth ) - len(dico[id_min_rigth ].strip('-'))
	seq_min_left = dico_aln_left[id_min_left]
	len_left = len(seq_min_left)
	dico_cut_rigth = {}
	dico_cut_left = {}
	for elt in dico.keys() :
		seq_left = dico[elt][0:len_left]
		if seq_left.count('-') < len(seq_left)*0.6 :
			dico_cut_left[elt] = seq_left
	for elt in dico.keys() :
		seq_rigth = dico[elt][len_rigth:len(dico[id_min_rigth])]
		if seq_rigth.count('-') < len(seq_rigth)*0.6 :
			dico_cut_rigth[elt] = dico[elt][len_rigth:len(dico[id_min_rigth])]

	with open(f"{prefix}_seq2.fasta",'w') as output_rigth, open(f"{prefix}_seq1.fasta",'w') as output_left:
		for elt in dico_cut_rigth.keys() : # Manque un control que ce soit bien un ATG !
			record = SeqRecord(Seq(dico_cut_rigth[elt].replace('-','')), id=str(elt), name=str(elt), description='')
			SeqIO.write(record,output_rigth, "fasta")
		for elt in dico_cut_left.keys():  # Manque un control que ce soit bien un ATG !
			record = SeqRecord(Seq(dico_cut_left[elt].replace('-','')), id=str(elt), name=str(elt), description='')
			SeqIO.write(record, output_left, "fasta")
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))