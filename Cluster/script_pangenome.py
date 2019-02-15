#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package script_pangenome.py
# @author Florian Charriat

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



########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens,comparaisonListe,isIn

if __name__ == "__main__":

	version="0.1"

############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used  to search orthologue groups of a sequence in fasta files''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display script_pangenome version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-g', '--group',type = str, required=True, dest = 'group', help = 'Path of the result of orthofinder (format csv)')
	filesreq.add_argument('-c', '--count',type = str, required=True, dest = 'count', help = 'Path of the count result of orthofinder (format csv)')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of the fasta which contains all sequence for pangenome analysis')
	filesreq.add_argument('--hote', type=str, required=True, dest='hote',
						  help='Path of the file containing species/host relationships')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdir', help = 'Path of the output directory')



######### Recuperation arguments ###########
	args = parser.parse_args()
	count = os.path.abspath(args.count)
	group= os.path.abspath(args.group)
	fasta = os.path.abspath(args.fasta)
	outdir =  os.path.abspath(args.outdir)
	hote=  os.path.abspath(args.hote)
	pathScript = sys.path[0]+'/Report_pangenome.Rmd'
########### Gestion directory ##############
	verifFichier(count)
	verifFichier(group)
	verifFichier(fasta)
	verifFichier(hote)
	outdir = verifDir(outdir)
	createDir([outdir])


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("      Welcome in secretome_Pipeline (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

	listeOG = []
	listeGene = list(fasta2dict(fasta).keys())
	with open(group,'r') as file_group :
		header = file_group.readline()
		for line in file_group :
			name = line.split()[0]
			genes = line.replace('\n','').replace('\t',',').split(':')[-1].split(',')[1:len(line.split())]
			if isIn(listeGene,genes) :
				listeOG.append(name)

	nb = 0
	with open(count,'r') as file_count,\
			open(f'{outdir}gene.count.csv','w') as result :
		header = f'Groupe{file_count.readline()}'.replace('\tTotal','')
		result.write(header)
		for line in file_count :
			name = line.split()[0]
			if name in listeOG :
				nb += 1
				total = line.split()[-1]
				gene = line.replace(f'\t{total}', '')
				result.write(gene)
				print(gene)
	# os.system("Rscript -e 'rmarkdown::render(%s, output_file=%s, quiet=TRUE, params = list(data = %s, dataH = %s))'"%('"'+pathScript+'"','"'+outdir+'Report.html'+'"','"'+outdir+'gene.count.csv"','"'+hote+'"'))

	############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))







