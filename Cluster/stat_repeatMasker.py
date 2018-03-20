#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The stat_repeatMasker script
	============================
	author: Charriat Florian\n
	contact: florian.charriat@inra.fr\n
	date: 13/03/2018\n
	version: 0.1

	Script description
	------------------

	his program is used to retrieve stat of all the data done by the script repeatMasker_build.py

	Example
	-------

	>>> stat_repeatMasker.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display stat_repeatMasker.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the result of the repeatMasker_build.py
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						Path of the output directory

"""


########## Module ###############
import argparse
import os

####### FUNCTION ################
def createDir(directory):
	''' Permet de vérifier si un dossier existe, si ce n'est pas le cas, 
	le dossier sera crée
	'''
	if not os.path.exists(directory):
		 	os.makedirs(directory)
	return

def verifDir(directory):
	'''
	Permet de mettre en forme le chemin du dossier pour être utilisé dans un script, 
	la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute
	'''
	if directory.endswith('/') == False :
		directory = directory +'/'

	return directory

if __name__ == "__main__":

########### Gestion directory #############
	directory = '/homedir/charriat/work/repeatMasker/result_repeatMasker/'
	outDir = '/homedir/charriat/work/repeatMasker/'
	directory = verifDir(directory)
	outDir = verifDir(outDir)

########## main script ######################
	listeID = []
	for ID in os.listdir(directory) :
		listeID.append(ID)
	listeID.sort()
	statAll = open(outDir+'stat_repeatMasker','w')
	statAll.write('Id souche'+'\t'+'Nombre séquence'+'\t'+'Longueur total'+'\t'+'Taux GC'+'\t'+'longueur bases masquée'+'\t'+ 'Pourcentage bases masquées'+'\n')
	statAll.close()
	for isolate in listeID :
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























