#!/usr/local/bioinfo/python/2.7.9/bin/python
# -*- coding: utf-8 -*-
# @package LDhat_launcher.py
# @author Florian CHARRIAT

"""
	The LDhat_launcher script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/07/2018
	:version: 0.1

	Script description
	------------------

	This program is used to generate stat of recombinaison with LDhat tools with fasta file input.
	
	Example
	-------

	>>> LDhat_launcher.py -f /homedir/user/work/file.fasta -o /homedir/user/work/result/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to//fasta/file/>, --fasta <path/to//fasta/file/>
						Path of the fasta to used
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
import argparse, egglib

pathLDhat = ''

##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	version = "0.1"

############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to generate stat of recombinaison with LDhat tools with fasta file input.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= 'display script_pangenome version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-dt', '--datatype', default=1, type=int, choices=[1, 2], dest='datatype',help='1 for haplotypic data (default), 2 for genotypic')
	filesreq.add_argument('-f', '--fasta', type=str, required=True, dest='fasta',
						  help='Path of the fasta to used')
	filesreq.add_argument('-o', '--outdir', type=str, required=True, dest='outdir', help='path of the output directory')

	######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta = os.path.abspath(args.fasta)
	outdir =  os.path.abspath(args.outdir)
	datatype = args.datatype
	outdir = verifDir(outdir)
	createDir([outdir])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print("\t" + form("|", 'yellow', 'bold') + form("        Welcome in LDhat_launcher (Version " + version + ")        ",type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	nameFasta = fasta.split('/')[-1].replace('.fasta','')

	##  Prétraitement du fichier fasta
	dico_fasta = fasta2dict(fasta)
	nbSeq = len(dico_fasta.keys())
	lenSeq = len(dico_fasta[dico_fasta.keys()[0]])

	with open('{}{}_LDhat.fasta'.format(outdir,nameFasta),'w') as new_fasta :
		new_fasta.write('{} {} {}\n'.format(nbSeq,lenSeq,datatype))
		for elt in dico_fasta.keys() :
			record = SeqRecord(seq=dico_fasta[elt].seq, id=elt, name=elt, description=dico_fasta[elt].description)
			SeqIO.write(record, new_fasta, "fasta")

	## Calcul du theta pour la création de la likehood table
	liste_seq = []
	for elt in dico_fasta.keys() :
		liste_seq.append([elt,str(dico_fasta[elt].seq),1])
	aln = egglib.Align.create(liste_seq)
	cs = egglib.stats.ComputeStats()
	cs.add_stats('lseff', 'thetaW')
	stats = cs.process_align(aln, max_missing=float(0.9))
	theta = round(stats['thetaW']/stats['lseff'],4)

	## Lancement LDhat
	print(form("\t- Launch of LDhat tools for convert the {} to locs ans sites file\n".format(fasta.split('/')[-1]),"green","bold"))
	os.system("{}/convert -seq {} -prefix {}{}_".format(pathLDhat,fasta,outdir,nameFasta))
	print(form("\t- Launch of LDhat tools for generate the Likehood table\n"))
	os.system("{}/complete -n {} -rhomax 100 -n_pts 101 -theta {} -prefix {}{}_".format(pathLDhat,nbSeq,theta,outdir,nameFasta))
	print(form("\t- Launch of LDhat tools for estimate the variable recombinaison rate\n","green","bold"))
	os.system("{}/interval -seq {}{}_sites.txt -loc {}{}_locs.txt -lk {}{}_new_lk.txt - prefix {}{} -its 10000000 -samp 5000 -bpen 5".format(pathLDhat,outdir,nameFasta,outdir,nameFasta,outdir,nameFasta,outdir,nameFasta))
	print(form("\t- Launch of LDhat tools for retrieve stat for the output of the last step\n", "green", "bold"))
	os.system("{}/stat -input {}{}_rates.txt -prefix {}{}_".format(pathLDhat, outdir,nameFasta,outdir,nameFasta))

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))