#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The Quality script
	==================

	author: CHARRIAT Florian\n
	contact: florian.charriat@inra.fr\n
	date: 9/03/2018\n
	version: 0.1

	Script description
	------------------

	This program is used to retrieve quality data of all the assemblies done by the script ABYSS_launch

	Example
	-------

	>>> Quality.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display ABYSS_launch.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the result of the ABYSS_launch.py
		- \-o <path/to/output/directory>, --outDirPath <path/to/output/directory>
						Path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import MODULES_SEB
from module_Flo import verifDir, createDir, form



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the fasta assembled files')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	outDir= args.outDirPath


if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display '+__file__+' version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the fasta files which must be assembled')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = args.dirPath
	outDir= args.outDirPath


########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)

	name_directory = [outDir]
	for folder in name_directory: 
		createDir(folder)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in ABYSS_launch (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################
	quality = open(outDir+"QualityAssembly.csv","w")
	quality.close()
	listeID = []
	qualityR = quality = open(outDir+"QualityR","w")
	qualityR.write("n\tn:500\tL50\tmin\tN80\tN50\tN20\tE-size\tmax\tsum\tname\t\n")
	qualityR.close()
	N50min = 25000
	exclus = open(outDir+"Exclus_Assembly","w")
	exclus.write("Liste des souches exclus, N50 <"+ str(N50min) +"\n\n")
	exclus.close()
	for ID in os.listdir(directory) :
		listeID.append(ID)
	listeID.sort()
	a = 0 
	nbAssemblage = 0
	for isolate in listeID : 	
		print("Recuperation des données Qualité de l'assemblage : " +isolate)
		if a == 0 :
			quality = open(outDir+"QualityAssembly.csv","a")
			quality.write("\nn\tn:500\tL50\tmin\tN80\tN50\tN20\tE-size\tmax\tsum\tname\t\n")
			quality.close()
		else : 
			Exclus = open(outDir+"Exclus_Assembly","a")
			Exclus.write("\nn\tn:500\tL50\tmin\tN80\tN50\tN20\tE-size\tmax\tsum\tname\t\n")
			Exclus.close()
		a = 0
		for kmers in [20,30,40,50,60,70,80,90]:
			nbAssemblage += 1
			for file in os.listdir(directory+isolate+"/abyss_assembly_"+isolate+"_"+str(kmers)):	
				if file.endswith('-stats.tab')==True :
					stat = open(directory+isolate+"/abyss_assembly_"+isolate+"_"+str(kmers)+"/"+file,"r")
					stats = stat.readlines()
					ligne = stats[3].split('\t')
					if int(ligne[6]) < N50min and a != 7: 
						Exclus = open(outDir+"Exclus_Assembly","a")
						Exclus.write(stats[3])
						Exclus.close()
						qualityR = open(outDir+"QualityR","a")
						qualityR.write(stats[3])
						qualityR.close()
						a = a + 1
						break
					if int(ligne[6]) < N50min and a == 7  :
						Exclus = open(outDir+"Exclus_Assembly","a")
						Exclus.write(stats[3])
						Exclus.close()
						qualityR = open(outDir+"QualityR","a")
						qualityR.write(stats[3])
						qualityR.close()
						break
					quality = open(outDir+"QualityAssembly.csv","a")
					quality.write(stats[3])
					qualityR = open(outDir+"QualityR","a")
					qualityR.write(stats[3])
					stat.close()
					quality.close()
					qualityR.close()



############## summary message #######################
	print(form('\n-------------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('- Quality a récupérées les données qualité des '+str(nbAssemblage)+' assemblages')
	print('- Les fichiers de sorties sont :\n')
	print('\t- QualityAssembly.csv : Fichier qui donne les statistiques des assemblages par souche')
	print('\t- Exclus_Assembly : Fichier qui donne les statistiques des assemblages exclus à cause de leur mauvaise qualité')
	print('\t- QualityR : Fichier contenant toutes les statistiques, ce fichier peut etre utilisé facilement avec R')
	print(form('-------------------------------------------------------------------------------------------------------------------------','red','bold'))


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




