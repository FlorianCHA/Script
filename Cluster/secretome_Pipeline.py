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

	This program is used to predict secretome with SignalP, TagetP and Phobius
	
	Example
	-------

	>>> ABYSS_launch.py -d /homedir/user/work/data -o /homedir/user/work/result

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
	filesreq.add_argument('-f', '--force', action='store_true', dest = 'force', help = 'Path of the output directory')

	
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

	for files in os.listdir(directory):
		if isFasta(files):
			fasta = open(directory+files,'r')
			lines = fasta.readlines()
			fasta.close()
			nb = 0
			part = 1
			outTargetP = '%s%s/%s_targetP.txt'%(outDir,recupId(files),recupId(files))
			outputPhobius = '%s%s/%s_phobius.txt'%(outDir,recupId(files),recupId(files))
			if os.path.exists('%s%s'%(outDir,recupId(files))) == True and force == False:
				raise ValueError(form("Le dossier output : '%s' existe deja, veuillez le supprimer ou utilisé la commande --force pour passer outre et supprimer les dossiers automatiquement"%(outDir+recupId(files)),"red","bold"))
			if os.path.exists('%s%s'%(outDir,recupId(files))) != 0 and force == True :
				os.system('rm -r %s%s'%(outDir,recupId(files)))
			createDir([outDir+recupId(files)+'/fasta_files',outDir+recupId(files)])
			os.system('touch %s %s' % (outTargetP,outputPhobius))
			listeTargetp = []
			listePhobius = []
	
			for line in lines: 
				if nb == 400 and line[0] == '>' :
					listePhobius.append('phobius.pl -short %s%s/fasta_files/%s_part%s.fasta  >> %s;\n'%(outDir,recupId(files),recupId(files),str(part),outputPhobius))
					listeTargetp.append('targetp -N %s%s/fasta_files/%s_part%s.fasta >> %s;\n'%(outDir,recupId(files),recupId(files),str(part),outTargetP))
					nb = 1
					part += 1
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,recupId(files),recupId(files),str(part)),'a')
					f.write(line)
					f.close
					
				
				elif line[0] == '>' :
					nb +=1
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,recupId(files),recupId(files),str(part)),'a')
					f.write(line)
					f.close
				
				else : 
					f = open('%s%s/fasta_files/%s_part%s.fasta'%(outDir,recupId(files),recupId(files),str(part)),'a')
					f.write(line)
					f.close
			
		
			f = open('%s/%s.sh'%(bash,recupId(files)),'w')
			f.write('#$ -o %s\n#$ -o %s\n#$ -N %s_secretome\n module load bioinfo/signalp/4.1\n'% (outDir+'error_files', outDir+'out_files',recupId(files)))
			f.write('\n########### Lancement targetP sur les fichiers fasta de 400 séquences ###################\n\n')
			for elt in listeTargetp :
				f.write(elt)
			f.write('\n########### Lancement Phobius sur les fichiers fasta de 400 séquences ###################\n\n')
			for elt in listePhobius :
				f.write(elt)
			f.close()







############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




