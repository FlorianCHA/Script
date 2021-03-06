#!/bin/env python
# -*- coding: utf-8 -*-
# @package Alignement.py
# @author Florian Charriat

"""
	The Alignement script
	=====================

	:author: CHARRIAT Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to align multiple RNA-seq on a genome and merge the different alignments into a file

	Example
	-------

	>>> Alignement.py -d /homedir/user/work/RNAseq -o /homedir/user/work/result -c /homedir/user/work/config.txt -r /homedir/user/work/RNAseq/assembly
	
	>>> Alignement.py -d /homedir/user/work/RNAseq -o /homedir/user/work/result -c /homedir/user/work/config.txt -r /homedir/user/work/RNAseq/assembly -f file.fasta
	

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Alignement.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the RNA-seq
		- \-r <path/to/reference/directory>, --refDir<path/to/reference/directory>
						path of directory that contains all the genome assembled
		- \-c <path/to/configFile>, --configFile <path/to/>configfile>
						path of the config file for toogle
		- \-o <path/to/output/directory>, --outDirPath <path/to/output/directory>
						Path of the output directory
	Input infos for running only one genome:
		- \-f <Name/file>, --file <Name/file>
						Input infos for running only one genome



"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFasta, isFastq ,  recupId



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to align multiple RNA-seq on a genome and merge the different alignments into a file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Alignement.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the RNA-seq')
	filesreq.add_argument('-r', '--refDir',type = str, required=True, dest = 'refDir', help = 'path of directory that contains all the genome assembled')
	filesreq.add_argument('-c', '--configFile',type = str, required=True, dest = 'configFile', help = 'Path of the config file for toogle')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')
	files = parser.add_argument_group('Input infos for running only one genome')
	files.add_argument('-f', '--file', type=str, required=False, dest = 'file', default ='', help = 'name file genome if the user only wants to process only one genome')
	
	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	ref = os.path.abspath(args.refDir)
	config = os.path.abspath(args.configFile)
	outDir= os.path.abspath(args.outDirPath)
	assembly = args.file

########### Gestion directory ##############
	directory = verifDir(directory,True)
	ref = verifDir(ref,True)
	verifFichier(config)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash/'
	braker = outDir+'bamForBraker/'
	name_directory = [outDir, bash,outDir+'error/',outDir+'output/',braker]
	for folder in name_directory: 
		createDir(folder)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("         Welcome in Alignement (Version " + version + ")           ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

########## Main ############################
	nbScript = 0
	run = open(outDir+'run_job_mapping.sh','w')
	run.close()
	for genome in os.listdir(ref) :
		if assembly in genome and isFasta(genome):
			nbScript += 1
			IDgenome = genome.replace('_scaffold','')
			IDgenome = recupId(IDgenome)
			print('\nCréation du script de mapping pour : ' + IDgenome)
			genomeOutDir = outDir+IDgenome
			resultMapping =	genomeOutDir+'/finalResults'
			nameFile = bash+IDgenome+'_mapping.sh'
			files = open(nameFile,'w')
			
			#Permet de géré les sortie de sge 
			files.write('#$ -o '+outDir+'output/'+IDgenome+'.out\n#$ -e '+outDir+'error/'+IDgenome+'.err\n')
			
			# Permet de charger puis lancer Toogle pour un alignement
			files.write('module load bioinfo/braker/1.9\nbioinfo/exonerate/2.4.7\nmodule load bioinfo/TOGGLE/0.3.6\n')
			files.write('rm -r '+genomeOutDir+'\n')
			files.write('toggleGenerator.pl -d '+directory+' -r '+ref+genome+' -c '+config+' -o '+genomeOutDir+' -nocheck;\n')
			# Permet de récupérer tous les hits accepté pour ensuite les merger
			files.write('cd '+resultMapping+';\nls *.accepted_hits.bam > bamList;\n')
			
			# Permet de merger les différents hits récupérés précédement
			mergefile = 'merged_'+IDgenome+'.accepted_hits.bam'
			files.write('samtools merge -b bamList -c '+mergefile+';\n')
			
			# Permet de triée les données du fichier bam contenant tous les mapping
			sortfile  = 'merged_'+IDgenome+'.accepted_hits_sort.bam'
			files.write('java -jar /usr/local/bioinfo/picard-tools/2.7.0//picard.jar SortSam I='+mergefile+' O='+sortfile+' SORT_ORDER=coordinate;\n')
			
			# Permet d'utiliser l'outils bam2hints pour formater les données pour l'annotation avec Augustus ou braker
			files.write('bam2hints --minintronlen=10 --maxintronlen=1000 --maxgaplen=9 --source=M --in='+sortfile+' --out=hints_'+IDgenome+'.raw.bam;\n')
			
			# Permet de selectionner seulement un set de read minimum requis pour un intron avec un script R
			files.write(sys.path[0]+'/filterHints.r -s '+IDgenome+' -p '+resultMapping+'\n')
			createDir(braker+IDgenome)
			files.write('cp -s '+resultMapping+'/hints_'+IDgenome+'.filtered.gff '+braker+IDgenome+'/')
			files.close()
			os.system('chmod 755 '+nameFile)
			run = open(outDir+'run_job_mapping.sh','a')
			run.write("qsub -N "+IDgenome+"_mapping -V -q long.q '"+nameFile+"'\n")
			run.close()



############## summary message #######################
	print(form('\n-------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\tInput :')
	print('\t\t- Donnée RNAseq : '+directory[:-1])
	print('\t\t- Repertoire génome : '+ref[:-1])
	if assembly != '' :
		print('\t\t- Génome : '+assembly)
	print('\t\t- Fichier config : '+config)
	
	print('\n\tOutput :')
	print('\t\t- script bash : '+bash[:-1])
	print('\t\t- fichier a lancer : '+outDir+'run_job_mapping.sh')	
	print('\t\t- Resultat des Jobs : '+braker[:-1])	
	
	print('\nSi vous souhaité lancer tous les '+str(nbScript)+' jobs de mapping veuillez taper la commande : ')
	print(form('\n\t\t\t\tbash '+outDir+'run_job_mapping.sh\n','green','bold'))
	print(form('-------------------------------------------------------------------------------------------------------------------','red','bold'))
	
		

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))














