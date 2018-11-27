#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Annotation_pipeline.py
# @author Florian Charriat

version = '0.2'
"""
	The Annotation_pipeline script
	==========================
	:author: CHARRIAT Florian
	:contact: florian.charriat@inra.fr
	:date: 12/11/2018
	:version: 0.2
	Script description
	------------------
	This program is used to annotate with BRAKER and augustus. This script use :\n\t
	- Toogle for align multiple RNA-seq on a genome and merge the different alignments into a file\n\t
	- Samtools for file is sort\n\t
	- Bam2hints for edit at hint format for braker\n\t
	- exonerate for align genome with annoted genome protein\n\t
	- Braker with HINTS option which use intron Hint of RNAseq and exonerate\n\t
	- Augustus with HINTS option which use intron and exon Hint of RNAseq and exonerate\n\t

	Example
	-------

	>>> Annotation_pipeline.py -d /homedir/user/work/RNAseq -o /homedir/user/work/result -c /homedir/user/work/config.txt -r /homedir/user/work/RNAseq/assembly p homedir/user/work/annotated_genome_protein.fasta -f file.fasta

	Help Programm
	-------------

	optional arguments:
	- \-h, --help
			show this help message and exit

	Input infos not mandatory:
	- \-v, --version
			Use if you want to know which version of Annotation_pipeline you are using
	- \-j NBJOB, --nbjob NBJOB
			Max number of job to run at the same time
	- \-s, --schema
			Use this option if you want a schema of the pipeline and not launch the Annotation_pipeline

	One of this input is mandatory:
	- \-r RNASEQ, --RNAseq RNASEQ
			Path of directory that contains all the rnaseq data.
			If you have only protein file, please use the -p
			option only.
	- \-p PROTEIN, --protein PROTEIN
			Path of annotated genome protein fasta file. If you
			have only RNAseq data, please use the -r option only.

	Input mandatory infos for running:
	- \-d DATADIR, --directory DATADIR
			Path of directory that contains all fasta file to annotate
	- \-S ID_SOUCHE, --specie ID_SOUCHE
			Name of the species (ex : magnaporthe_oryzae)
	- \-o OUTDIRPATH, --outDir OUTDIRPATH
			Path of the output directory
"""

########## Module ###############
## Python modules
import argparse, os, sys, re
from argparse import RawTextHelpFormatter

# Import module_Flo
from script.module_Flo import verifDir, createDir, form, verifFichier, isFasta, isFastq, recupId


descriptionTools = '''This program is used to annotate with BRAKER and augustus. This script use :\n\t
	- Toogle for align multiple RNA-seq on a genome and merge the different alignments into a file\t
	- Samtools for file is sort\t
	- Bam2hints for edit at hint format for braker\t
	- exonerate for align genome with annoted genome protein\t
	- Braker with HINTS option which use intron Hint of RNAseq and exonerate\t
	- Augustus with HINTS option which use intron and exon Hint of RNAseq and exonerate\t'''

epilogTools = ""

def mandatoryDatas(x):
	if parserOther.parse_known_args()[0].RNAseq == None and parserOther.parse_known_args()[0].Protein == None :
		raise parserOther.error(form('\n\n\tError : One of the following arguments are required: -r/--RNASeq and/or -p/--protein\n\tYou must give at least a protein file or a RNAseq data to proceed the annotation\n','red','bold'))
	return x

class InputDir(object):
	"""For check initial dir"""

	def __init__(self, path = None):
		"""Initialise variable"""
		#Change relative path to absolute path
		self.path = os.path.abspath(path)
		self.directoryList = []
		self.filesList = []
		if self.path[-1] != "/":
			self.path += "/"
		if not os.path.isdir(self.path):
			raise argparse.ArgumentTypeError("\n\n\tERROR : path '{}' is not valid path.\n\n\n\nExiting...".format(self.path ))
		if not os.access(self.path, os.W_OK):
			raise argparse.ArgumentTypeError("\n\n\tERROR : Permission denied for our path '{}' check admin system.\n\n\n\nExiting...".format(self.path ))

		# appel les fonctions
		self.testDirContent()
		self.testAllowExtention()

	def testDirContent(self):
		"""List all in directory"""
		self.listPath = [entry.name for entry in os.scandir(self.path)]
		for fileIn in os.scandir(self.path):
			if fileIn.is_dir():
				self.directoryList.append(fileIn.name)
			if fileIn.is_file():
				self.filesList.append(fileIn.path)

	def testAllowExtention(self):
		""""""
		findExtentionList = []
		for fileIN in self.filesList:
			## file are fasta
			if re.findall("\.fasta$|fa$|fna$|faa$",fileIN, flags=re.IGNORECASE):
				extension = fileIN.split(".")[-1]
				findExtentionList.append(extension)
				if extension not in ["fasta"]:
					goodName = fileIN.replace(extension, "fasta")
					os.rename(fileIN, goodName)
					print(form("\n\tWarnings :Rename {} to {}\n".format(fileIN, goodName), 'yellow','bold'))

		if len(findExtentionList) == 0 :
			raise argparse.ArgumentTypeError('\n\n\tERROR :no fasta, fa or fna file found .\n\n\nExiting...')

	def __repr__(self):
		return "{}\n\n{}".format(self.__class__, self.__dict__)

	def __str__(self):
		"""Fonction qui permet de formater le text de sortie lors du print du dictionnaire"""
		return """path = \t\t{0}\n\nfiles/dir = \t{1}""".format(self.path,  "\n\t\t".join(self.listPath))


if __name__ == "__main__":

############ Argparse #####################
	parserOther = argparse.ArgumentParser(prog=__file__,
									 formatter_class=argparse.RawDescriptionHelpFormatter,
									 description=descriptionTools,
									 epilog=epilogTools)

	inOutOptional = parserOther.add_argument_group('Input infos not mandatory')
	inOutOptional.add_argument('-v', '--version', action='version', version=version, help='Use if you want to know which version of Annotation_pipeline you are using')
	inOutOptional.add_argument('-j', '--nbjob', type=int, required=False, default=100, dest='nbJob', help='Max number of job to run at the same time')
	inOutOptional.add_argument('-s', '--schema',action='store_true', dest='schema', help='Use this option if you want a schema of the pipeline and not launch the Annotation_pipeline')
	inOutOptional.add_argument('-R', '--autorun',action='store_true', dest='autorun', help='Use this option if you want a run the Annotation_pipeline')

	inOutMandatory2 = parserOther.add_argument_group('One of this input is mandatory')
	inOutMandatory2.add_argument('-r', '--RNAseq', type=str, required=False,default=None, dest='RNAseq', help='Path of directory that contains all the rnaseq data. If you have only protein file, please use the -p option only.')
	inOutMandatory2.add_argument('-p', '--protein', type=str, required=False,default=None,dest='Protein', help='Path of annotated genome protein fasta file. If you have only RNAseq data, please use the -r option only.')

	parserMandatory = argparse.ArgumentParser(
									 parents=[parserOther],
									 add_help=False
									)
	inOutMandatory = parserMandatory.add_argument_group('Input mandatory infos for running')
	inOutMandatory.add_argument('-d', '--directory', type=InputDir, default=None, required=True, dest='datadir', help='Path of directory that contains all fasta file to annotate')
	inOutMandatory.add_argument('-S', '--specie', type=str, required=True, dest='id_souche', help='Name of the species (ex : magnaporthe_oryzae)')
	inOutMandatory.add_argument('-o', '--outDir', type=mandatoryDatas, required=True, dest='outDirPath', help='Path of the output directory')

############ Recupération Argument #####################

	args = parserMandatory.parse_args()

	# print(args)
	datadir = args.datadir.path
	RNAseq = args.RNAseq
	nbJob  = args.nbJob
	schema = args.schema
	Protein = args.Protein
	id_souche = args.id_souche
	outDir= args.outDirPath

	################################################################################
	#Control des inputs
	outDir = verifDir(outDir)
	createDir([outDir])
	if RNAseq != None :
		RNAseq = os.path.abspath(RNAseq)
		RNAseq = verifDir(RNAseq,check = True)
	else :
		print(form("\n\tWarnings : you didn't give a directory for RNAseq data, the annotation will be done only with the protein file\n", 'yellow','bold'))

	if Protein != None :
		Protein = os.path.abspath(Protein)
		verifFichier(Protein)
	else :
		print(form("\n\tWarnings : you didn't give a protein file, the annotation will be done only with the directory for RNAseq data\n", 'yellow','bold'))

	id_souche = id_souche.replace(' ','_')

	confTXT = '''
{{
		"datadir"	: "{}",
		"OutDir"	: "{}",
		"RnaSeqDir"	: "{}",
		"id_souche"	: "{}",
		"protRef"	: "{}"
}}
	'''.format(datadir, outDir, RNAseq, id_souche, Protein)

	with open('./SupplementaryFile/config.json','w') as f:
		f.write(confTXT) # Ecriture de l'entete du fichier de configuration

	if schema :
		with open('makeFig.sh','w') as f:
			f.write('#!/bin/bash\n')
			f.write('module unload system/python/2.7.9\n')
			f.write('module load bioinfo/snakemake/3.13.3\n')
			f.write('snakemake -s BRAKER_pipeline2.snake --dag | dot -Tpdf > schema_pipeline.pdf')
		print(form('\nCreating the pipeline schema in schema_pipeline.pdf','green','bold'))
		print(form('-'*len('Creating the pipeline schema in schema_pipeline.pdf')+'\n','yellow','bold'))
		os.system('bash makeFig.sh')
		os.system('rm makeFig.sh')

	else :
		with open('Launcher.sh', 'w') as f:
			f.write('#!/bin/bash\n')
			f.write('module purge\n')
			f.write('module load bioinfo/snakemake/3.13.3\n')
			f.write('snakemake -s BRAKER_pipeline2.snake --jobs {0} --cluster "qsub -q long.q -cwd -V -pe parallel_smp {{threads}} -l mem_free={{params.l_mem_free}}"'.format(nbJob))

		print(form('\nLaunch of Annotation pipeline','green','bold'))
		print(form('-'*len('Launch of Annotation pipeline'),'yellow','bold'))
		print(form('qsub -V -cwd -N ANNOT -q long.q -b Y ./Launcher.sh\n\n','white','bold'))

		if args.autorun:
			os.system('qsub -V -cwd -N ANNOT -q long.q -b Y ./Launcher.sh')
		# os.system('rm Launcher.sh')
