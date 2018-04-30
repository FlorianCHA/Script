#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package secretome_Pipeline.py
# @author Florian Charriat

"""
	The secretome_Pipeline script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/04/2018
	:version: 0.1

	Script description
	------------------

	This program is used to predict secretome with SignalP, TagetP, Phobius.\n
	This program uses the comparaisonSecretome script to retrieve and compare information from the secretome prediction tools.\n
	This program uses also the selection_TMHMM script to select, from TMHMM ouput, only protein with no TH domain or only one TH in 60 first aa.\n
	Please make sure that this script is present in the directory.
	
	Example
	-------

	>>> secretome_Pipeline.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fasta files
						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory
		- \-f                             --force
						force the script to remove the old output data
"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to predict secretome with SignalP, TagetP and Phobius''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display secretome_Pipeline version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the fasta files')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')
	filesreq.add_argument('-f', '--force', action='store_true', dest = 'force', help = 'force the script to remove output data')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outDir= os.path.abspath(args.outdirPath)
	force = args.force

########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash'
	name_directory = [outDir,outDir+'error_files', outDir+'out_files',bash]
	createDir(name_directory)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in ABYSS_launch (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################
	nbfile = 0
	for files in os.listdir(directory):
		if isFasta(files):
			idFile = recupId(files)
			nbfile +=1
			fasta = open(directory+files,'r')
			lines = fasta.readlines()
			fasta.close()
			nb = 0 # Permet d'initialiser une variable qui servira a séparer le fichier fasta en plusieurs fichier fasta
			nbSeq = 0 # Permet de savoir le nombre de séquence presente dans le fichier fasta
			part = 1 # Permet de nomer les fichiers avec un numero part morceau de sequence
			
			########## Création dossier de résultat ############
			if os.path.exists('%s%s'%(outDir,idFile)) == True and force == False:
				raise ValueError(form("Le dossier output : '%s' existe deja, veuillez le supprimer ou utilisé la commande --force pour passer outre et supprimer les dossiers automatiquement"%(outDir+idFile),"red","bold"))
			if os.path.exists('%s%s'%(outDir,idFile)) != 0 and force == True :
				os.system('rm -r %s%s'%(outDir,idFile))
			outDirFasta = outDir+idFile+'/fasta_files'
			outDir_comparaison = outDir+idFile+'/1_predicted/'
			outDir_selectTMHMM = outDir+idFile+'/2_selectTMHMM/'
			name_directory = [outDir+idFile,outDirFasta,outDir_comparaison,outDir_selectTMHMM]
			createDir(name_directory)
			
			######### Création de fichier de résultat ##########
			outTargetP = '%s%s_targetP.txt'%(outDir_comparaison,idFile)
			outputPhobius = '%s%s_phobius.txt'%(outDir_comparaison,idFile)
			outputSignalP = '%s%s_signalP.txt'%(outDir_comparaison,idFile)
			os.system('touch %s %s %s' % (outTargetP,outputPhobius,outputSignalP))

			listeTargetp = []
			listePhobius = []
			listeSignalP = []
	
			for line in lines: 
				if nb == 400 and line[0] == '>' :
					listePhobius.append('phobius.pl -short %s%s/fasta_files/%s_part%s.fasta  >> %s;\n'%(outDir,idFile,idFile,str(part),outputPhobius))
					listeTargetp.append('targetp -N %s%s/fasta_files/%s_part%s.fasta >> %s;\n'%(outDir,idFile,idFile,str(part),outTargetP))
					listeSignalP.append('signalp -u 0.34 -U 0.34 %s%s/fasta_files/%s_part%s.fasta >> %s;\n'%(outDir,idFile,idFile,str(part),outputSignalP))
					nb = 1
					nbSeq +=1
					part +=1
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,idFile,idFile,str(part)),'a')
					f.write(line)
					f.close
					
				
				elif line[0] == '>' :
					nb +=1
					nbSeq +=1
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,idFile,idFile,str(part)),'a')
					f.write(line)
					f.close
				
				else : 
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,idFile,idFile,str(part)),'a')
					f.write(line)
					f.close
					outputSignalP
					
			# Les deux prechaines lignes permet de traiter le dernier fichier qui fait moint de 400 sequences
			if nbSeq%400 != 0 :
				listePhobius.append('phobius.pl -short %s%s/fasta_files/%s_part%s.fasta  >> %s;\n'%(outDir,idFile,idFile,str(part),outputPhobius))
				listeTargetp.append('targetp -N %s%s/fasta_files/%s_part%s.fasta >> %s;\n'%(outDir,idFile,idFile,str(part),outTargetP))
				listeSignalP.append('signalp -u 0.34 -U 0.34 %s%s/fasta_files/%s_part%s.fasta >> %s;\n'%(outDir,idFile,idFile,str(part),outputSignalP))
					
			
			f = open('%s/%s_secretomeTools.sh'%(bash,idFile),'w')

			f.write('\n########### Lancement listeSignalP sur les fichiers fasta de 400 séquences ###################\n\n')
			for elt in listeSignalP :
				f.write(elt)
			f.write('\n########### Lancement targetP sur les fichiers fasta de 400 séquences ###################\n\n')
			for elt in listeTargetp :
				f.write(elt)
			f.write('\n########### Lancement Phobius sur les fichiers fasta de 400 séquences ###################\n\n')
			for elt in listePhobius :
				f.write(elt)
			f.close()
			f = open('%s/%s.sh'%(bash,idFile),'w')
			f.write('#$ -e %s\n#$ -o %s\n#$ -N %s_secretome\n#$ -q normal.q\n#$ -V\n\nmodule load bioinfo/signalp/4.1\n\n'% (outDir+'error_files', outDir+'out_files',idFile))
			f.write('\n\n%s Lancement des trois outils de prédiction %s\n\n'%("#"*10,"#"*10))
			f.write('bash %s/%s_secretomeTools.sh\n\n'%(bash,idFile))
			f.write('%s Comparaison des 3 outils de prédiction %s\n\n'%("#"*10,"#"*10))
			f.write('comparaisonSecretome.py -o %s --phobius %s --targetp %s --signalp %s --rank 2 --fasta %s%s\n'%(outDir_comparaison,outputPhobius,outTargetP,outputSignalP,directory,files))
			f.write('\n\n%s Lancement TMHMM %s\n\n'%("#"*10,"#"*10))
			f.write('tmhmm -short %s%s_secreted_1.fasta > %s%s_TMHMM.txt\n'%(outDir_comparaison,idFile,outDir_selectTMHMM,idFile))
			f.write('\n\n%s Selection des proteines en fonction du TMHMM %s\n\n'%("#"*10,"#"*10))
			f.write('selection_TMHMM.py -t %s%s_TMHMM.txt -f  %s%s_secreted_1.fasta -o %s\n'%(outDir_selectTMHMM,idFile,outDir_comparaison,idFile,outDir_selectTMHMM))
			print(form('Script créé pour %s\n'%idFile,'green','bold'))
############## summary message #######################

	print(form('\n-----------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\n\tInput : \n\t\t- '+ directory[:-1])
	print('\n\tOutput :')
	print('\t\t - script bash créé : ' +bash)
	print('\t\t - Résultat des prédictions des secretomes : '+outDir)
	print('\nThe secretome_Pipeline a traité les '+str(nbfile)+' fichiers fasta, veuillez taper la commande suivante pour lancer les scripts crés :\n\n\t\t\t\t'+form('qsub %s/%s.sh\n'%(bash,idFile),'green','bold'))
	print(form('----------------------------------------------------------------------------------------------------------------------','red','bold'))





############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




