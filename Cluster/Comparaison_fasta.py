#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Comparaison_fasta.py
# @author Florian CHARRIAT

"""
	The Comparaison_fasta script
	===========================
	:author: Florian CHARRIAT
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to compare compare the content of two fasta (name of sequences or sequences).
	Example
	-------
	>>> Comparaison_fasta.py -f1 /homedir/file1.fasta -f2 /homedir/file2.fasta -o result/ --name


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Comparaison_fasta.py version number and exit
	Input mandatory infos for running:
		- \--file1 <path/to/fasta/file>
						path of the first fasta file to used
		- \--file2 <path/to/fasta/file>
						path of the second fasta file to used
		- \-o <path/to/output/directory/>, --output <path/to/output/directory/>
						path of the output directory
	Input infos for running with default values:
		- \--database1 <path/to/fasta/file/>
						If file1 are alignement file, use this option and give the fasta file used for alignement
		- \--database2 <path/to/fasta/file/>
						If file2 are alignement file, use this option and give the fasta file used for alignement
		- \-n , --name
						compare name of sequence between fasta file
		- \-seq , --sequence
						compare name of sequence between fasta file (use blast for the sequence compraison)
"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human, comparaisonListe

## Python modules
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	version = '0.1'
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to compare compare the content of two fasta (name of sequences or sequences).''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help='display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument( '--file1',type = str, required=True, dest = 'fasta1', help = 'Path of the first fasta file to used')
	filesreq.add_argument( '--file2',type = str, required=True, dest = 'fasta2', help = 'Path of the second fasta file to used')


	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-n', '--name',action='store_true',  dest = 'name',help = 'Compare name of sequence between fasta file')
	files.add_argument('-o', '--output', type=str,  required=False, default = 'None', help='Path of the output directory (option used only for sequence comparison)')
	files.add_argument('-seq', '--sequence',action='store_true', dest = 'seq', help = 'Compare name of sequence between fasta file (use blast for the sequence comparison)')
	filesreq.add_argument('--database1', type=str, required=False, default = 'None',dest='db1',
						  help='If file1 is alignement file, use this option and give the fasta file used for alignement (for seq option only)')
	filesreq.add_argument('--database2', type=str, required=False,  default = 'None', dest='db2',
						  help='If file2 is alignement file, use this option and give the fasta file used for alignement (for seq option only)')
######### Recuperation arguments ###########
	args = parser.parse_args()
	if args.output != 'None' :
		output = os.path.abspath(args.output)
	else :
		output =args.output
	fasta1 = os.path.abspath(args.fasta1)
	fasta2 = os.path.abspath(args.fasta2)
	name = args.name
	seq = args.seq
	db1 = args.db1
	db2 = args.db2

########### Gestion directory ##############
	if name == False and seq == False :
		print(form('\nError : Any option of comparison given, please choose the option seq (--sequence) and/or name (--name)\n','red','bold'))
		exit()
	if seq == True and output == 'None' :
		print(form("\nError : You didn't give any output while you chose the comparison, please use the --output option for give a output directory\n",'red', 'bold'))
		exit()
	verifFichier(fasta1)
	verifFichier(fasta2)
	if db1 != 'None' :
		verifFichier(db1)
		db1 = os.path.abspath(db1)
	if db2 != 'None' :
		verifFichier(db2)
		db2 = os.path.abspath(db2)

	output =  verifDir(output)
	listeDir = [output]


	if db1 != 'None' or db2 != 'None' :
		dfasta = output + 'fasta_file_from_aln/'
		listeDir.append(dfasta)
	if seq :
		Bresult = output+'Blast/'
		listeDir.append(Bresult)
		Bresult1 = output + 'Blast/1_blast/'
		listeDir.append(Bresult1)


	createDir(listeDir)

############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("               Comparaison_fasta (Version " + version + ")         ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################

	liste1 = list(fasta2dict(fasta1).keys())
	liste2 = list(fasta2dict(fasta2).keys())
	if len(liste2)< len(liste1) :
		fasta_inter = fasta1
		liste_inter = liste1
		db_inter = db1
		fasta1 = fasta2
		liste1 = liste2
		db1 = db2
		fasta2 = fasta_inter
		liste2 = liste1
		db2 = db_inter


	nameF1 = fasta1.split('/')[-1]
	nameF2 = fasta2.split('/')[-1]

	if name:
		list = comparaisonListe(liste1, liste2)
		print(form(f'They have {len(list)} sequences name in common between {fasta1.split("/")[-1 ]} (nb sequence : {len(liste1)}) and {fasta2.split("/")[-1 ]} (nb sequence : {len(liste2)})\n','green','bold'))


	if seq :
		if db1 !='None' :
			print(form(f'Warning : you use the database1 option, the pipeline recovers the sequences from the file : {db1}\n','orange','bold'))
			with open(f'{dfasta}{nameF1}.fasta','w') as f :
				fasta = fasta2dict(db1)
				aln = fasta2dict(fasta1)
				for elt in sorted(aln.keys(), key = sort_human) :qsta
					sequence = fasta[elt].seq
					record = SeqRecord(sequence, id=str(elt), name=str(elt), description=fasta[elt].description)
					SeqIO.write(record, f, "fasta")
			fasta1 = f'{dfasta}{nameF1}.fasta'
		if db2 !='None' :
			print(form(f'Warning : you use the database2 option, the pipeline recovers the sequences from the file : {db2}\n','orange','bold'))
			with open(f'{dfasta}{nameF2}.fasta','w') as f :
				fasta = fasta2dict(db2)
				aln = fasta2dict(fasta2)
				for elt in sorted(aln.keys(), key = sort_human) :
					sequence = fasta[elt].seq
					record = SeqRecord(sequence, id=str(elt), name=str(elt), description=fasta[elt].description)
					SeqIO.write(record, f, "fasta")
			fasta2 = f'{dfasta}{nameF2}.fasta'

		db_fasta2 = f"{Bresult}0_database/{nameF2.replace('.fasta','')}"
		output_blast1= f"{Bresult}1_blast/{nameF2.replace('.fasta','')}_blast_result.txt"
		print(form('\t- Création de la base de données pour psiblast','white','bold'))
		os.system(f'makeblastdb -dbtype prot -in {fasta2} -input_type fasta -out {db_fasta2}\n')
		print(form(f'\n\t- Lancement de blast (query : {nameF1} , db : {nameF2.replace(".fasta","")})','white','bold'))
		os.system(f'\nblastp -query {fasta1} -db {db_fasta2} -evalue 1e-3 -out {output_blast1}')
		print(form(f'\t- Utilisation du fichier de sortie de blast pour la comparaison\n','white','bold'))
		listeID = []
		listeFalse = []

		with open(output_blast1,'r') as f :
			for line in f :
				if line[0:6] == 'Query=' :
					id = line.split()[1]
				elif line[0:7] == 'Length=' :
					length = int(line.replace('\n','').split('Length=')[1])
				elif 'Identities' in line  :
					prcIdentity = int(line.split('(')[1].split('%')[0])
					couverture = int(line.split('/')[-1].split()[0])
					if couverture > length*0.8 and prcIdentity > 80 :
						if id not in listeID :
							listeID.append(id)
					else :
						if id not in listeFalse and id not in listeID :
							listeFalse.append(id)

		print(form(f'They have {len(listeID)} sequences in common between {nameF1} (nb sequence : {len(liste1)}) and {nameF2} (nb sequence : {len(liste2)})\n','green','bold'))

		dico_file1 = fasta2dict(fasta1)
		with open(f'{output}notInF2.fasta','w') as f :
			for elt in sorted(listeFalse,key = sort_human) :
				sequence = dico_file1[elt].seq
				record = SeqRecord(sequence, id=str(elt), name=str(elt), description=dico_file1[elt].description)
				SeqIO.write(record, f, "fasta")

		with open(f'{output}common_sequence.fasta','w') as f :
			for elt in sorted(listeID,key = sort_human) :
				sequence = dico_file1[elt].seq
				record = SeqRecord(sequence, id=str(elt), name=str(elt), description=dico_file1[elt].description)
				SeqIO.write(record, f, "fasta")



############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))