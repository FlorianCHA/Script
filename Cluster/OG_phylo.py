#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package OG_phylo.py
# @author Florian CHARRIAT

"""
	The script_pangenome script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/07/2018
	:version: 0.1

	Script description
	------------------

	This program is used to search orthologue groups of a sequence in fasta files
	
	Example
	-------

	>>> script_pangenome.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-a <path/to/alingment/HMM/file>, --align <path/to/fasta/file>
						path of the alignement HMM file of HmmerPipeline.py
		- \-db <path/to/blast/database>, --database <path/to/prosite.dat/file>
						path of the database of all genome which must be contain the database crate with makeblastbd.
		- \-f <path/to/secretome/fasta/file/>, --fasta <path/to/prosite.dat/file>
						path of the fasta file of all secretome
		- \-g                             --orthoGroups
						path of the directory which contain OrthoFinder result.
						
		- \-p                             --phenotype
						path of the phenotype file. Please, change the use the model for create your phenotype file.						
					
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
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to search orthologue groups of core genome single copy and create a consensus for phylogenetic tree''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= 'display script_pangenome version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-g', '--group', type=str, required=True, dest='group',
						  help='Path of the result of orthofinder (format csv)')
	filesreq.add_argument('-c', '--count', type=str, required=True, dest='count',
						  help='Path of the count result of orthofinder (format csv)')
	filesreq.add_argument('-f', '--fasta', type=str, required=True, dest='fasta',
						  help='Path of the fasta which contains all sequence used for orthofinder')
	filesreq.add_argument('-o', '--outdir', type=str, required=True, dest='outdir', help='Path of the output directory')

	######### Recuperation arguments ###########
	args = parser.parse_args()
	count = os.path.abspath(args.count)
	group= os.path.abspath(args.group)
	fasta = os.path.abspath(args.fasta)
	outdir =  os.path.abspath(args.outdir)
	outdir = verifDir(outdir)
	tranDir = outdir+'translatorX_output/'
	createDir([outdir,tranDir])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print("\t" + form("|", 'yellow', 'bold') + form("            Welcome in OG_phylo (Version " + version + ")         ",type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	nb = 0
	listeOG = []
	print(form(f'\t- Traitement du fichier {count} pour récupérer les groupes orthologue single copy core genome','white','bold'))
	print()
	# with open(count, 'r') as file_count:
	# 	header = file_count.readline()
	# 	for line in file_count :
	# 		name = line.split()[0]
	# 		nbGenes = line.split()[1:len(line.split())]
	# 		if nbGenes.count('1') == (len(nbGenes) -1 ) :
	# 			nb += 1
	# 			listeOG.append(name)

	print(form(f'\t- {nb} groupe orthologue single copy core genome on été trouvé','white','bold'))
	print()
	print(form(f'\t- Traitement du fichier {group} pour récupérer les sequences','white','bold'))

	# dico = fasta2dict(fasta)
	# with open(group, 'r') as file_group:
	# 	header = file_group.readline()
	# 	for line in file_group :
	# 		name = line.split()[0]
	# 		if name in listeOG :
	# 			Genes = line.split()[1:len(line.split())]
	# 			with open(f'{outdir}{name}.fasta','w') as f :
	# 				for gene in Genes :
	# 					sequence = dico[gene].seq
	# 					record = SeqRecord(sequence, id=str(gene), name=str(gene),description= dico[gene].description)
	# 					SeqIO.write(record, f, "fasta")
	print()
	print(form(f'\t- Lancement de transletorX', 'white', 'bold'))
	# for elt in os.listdir(outdir) :
	# 	if elt.endswith('.fasta') :
	# 		print(f'{elt.replace(".fasta","")} in process')
	# 		os.system(f'translatorx_vLocal.pl -i {outdir}{elt} -o {tranDir}{elt.replace(".fasta","")}> stdout.txt')

	dico_seq = {}
	for file_name in os.listdir(tranDir):
		if file_name.endswith('.nt_ali.fasta') :
			file_path = f'{tranDir}{file_name}'
			dico = fasta2dict(file_path)
			dico_species = {elt.split('_')[1] : value for elt,value in dico.items()}
			for elt in sorted(dico_species.keys(), key = sort_human) :
				if elt not in dico_seq.keys() :
					dico_seq[elt] = str(dico_species[elt].seq)
				else :
					dico_seq[elt] = dico_seq[elt] + str(dico_species[elt].seq)

	print()
	print(form(f'\t- Merged of all fasta', 'white', 'bold'))
	with open(f'{outdir}sequence_merged.fasta','w') as f :
		for elt in sorted(dico_seq.keys(), key = sort_human) :
			sequence = Seq(dico_seq[elt])
			record = SeqRecord(sequence, id=str(elt), name=str(elt),description= f'| length : {len(sequence)}')
			SeqIO.write(record, f, "fasta")
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))