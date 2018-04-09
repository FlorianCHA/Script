#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The stat_repeatMasker script
	============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 13/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to retrieve the statistique of result from repeatMasker_build.py script

	Example
	-------

	>>> stat_repeatMasker.py -d /homedir/user/work/repeatMasker_result/ -o /homedir/user/work/stat_repeatMasker

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display stat_repeatMasker.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the result of the repeatMasker_build.py (output + repeatMasker_result)
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						Path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier


if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''TThis program is used to retrieve the statistique of result from repeatMasker_build.py script''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display stat_repeatMasker.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the result of the repeatMasker_build.py (output + repeatMasker_result')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outDir= os.path.abspath(args.outdirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	name_directory = [outDir]
	createDir(name_directory)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("      Welcome in stat_repeatMasker (Version " + version + ")       ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')
	
	
	
########## main script ######################
	listeID = []
	for ID in os.listdir(directory) :
		listeID.append(ID)
	listeID.sort()
	nbAssemblage = 0
	statAll = open(outDir+'stat_repeatMasker.txt','w')
	statAll.write('Id souche'+'\t'+'Nombre séquence'+'\t'+'Longueur total'+'\t'+'Taux GC'+'\t'+'longueur bases masquée'+'\t'+ 'Pourcentage bases masquées'+'\n')
	statAll.close()
	for isolate in listeID :
		nbAssemblage += 1
		print(form("Récupération des statistique de l'assemblage : "+isolate,'green'))
		for files in os.listdir(directory+isolate):	
			if files.endswith('.fa.tbl')==True :
				filePath = directory+isolate+'/'+files
				stat = open(filePath,'r')
				stat = stat.readlines()
				Nbseq = stat[2].split(':')
				total_length = stat[3].split(':')
				total_length = total_length[1].split('(')
				GC_level = stat[4].split(':')
				bases_masked = stat[5].split(':')
				bases_masked = bases_masked[1].split('(')
				statAll = open(outDir+'stat_repeatMasker','a')
				statAll.write(files.replace('.fa.tbl','')+'\t'+ Nbseq[1].strip()+'\t'+total_length[0].strip()+'\t'+GC_level[1].strip()+'\t'+bases_masked[0].strip() +'\t'+bases_masked[1].replace(')','').strip()+'\n')
				statAll.close()




############## summary message #######################

	print(form('\n-------------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))

	print('\tInput :')
	print('\t\t- Repertoire des resultats de repeatMasker_build : '+directory[:-1])	
	
	print('\n\tOutput :')
	print('\t\t- Résultat : '+outDir[:-1])
	
	print('\n\tFichier de sortie :')
	print("\t\t- stat_repeatMasker.txt : Fichier qui donne les statistiques de l'outils repeatMasker des"+str(nbAssemblage)+' assemblages\n')
	
	print('\n Quality a récupérées les données qualité des '+str(nbAssemblage)+' assemblages\n')

	print(form('-------------------------------------------------------------------------------------------------------------------------','red','bold'))



############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))



















