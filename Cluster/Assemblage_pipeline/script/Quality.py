#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The Quality script
	==================

	:author: CHARRIAT Florian\n
	:contact: florian.charriat@inra.fr\n
	:date: 9/03/2018\n
	:version: 0.1

	Script description
	------------------

	This program is used to retrieve quality data of all the assemblies done by the script ABYSS_launch

	Example
	-------

	>>> Quality.py -d /homedir/user/work/data/result/ -o /homedir/user/work/quality

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Quality.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta assembled by ABYSS_launch (output + result)
		- \-o <path/to/output/directory>, --outFilePath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to retrieve the quality data of the assembly done by the ABYSS_launch script''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Quality.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the fasta assembled by ABYSS_launch (output + result)')
	filesreq.add_argument('-o', '--outFile',type = str, required=True, dest = 'outFilePath', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outFile= os.path.abspath(args.outFilePath)


if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to retrieve the quality data of the assembly done by the ABYSS_launch script''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the assembly files which assembled')
	filesreq.add_argument('-o', '--outFile',type = str, required=True, dest = 'outFilePath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	outFile= args.outFilePath


########### Gestion directory ##############
	directory = verifDir(directory,True)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in    Quality   (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################
	listeID = []
	for ID in os.listdir(directory) :
		listeID.append(ID)
	listeID.sort()
	
	qualityR = quality = open(outFile,"w")
	qualityR.write("n\tn:500\tL50\tmin\tN80\tN50\tN20\tE-size\tmax\tsum\tname\t\n")
	for isolate in listeID : 	
		for kmers in [20,30,40,50,60,70,80,90]:
			for file in os.listdir(directory+isolate+"/"+isolate+"_"+str(kmers)):	
				if file.endswith('-stats.tab')==True :
					stat = open(directory+isolate+"/"+isolate+"_"+str(kmers)+"/"+file,"r")
					stats = stat.readlines()
					ligne = stats[3].split('\t')
					qualityR.write(stats[3])

	qualityR.close()
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




