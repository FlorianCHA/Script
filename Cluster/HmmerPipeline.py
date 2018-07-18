#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package HmmerPiepline.py
# @author Florian Charriat

"""
	The HmmerPiepline script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 27/06/2018
	:version: 0.1

	Script description
	------------------

	This program is used to search a HMM profil for search EffecteurMAX. This program used HMMbuild et HMMsearch.

	Example
	-------

	>>> HmmerPiepline.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display HmmerPiepline.py version number and exit
						
	Input mandatory infos for running:
		- \-a <path/to/alignement/file>, --alignementFile <path/to/alignement/file>
						path of alignement file contrain by structure
		- \-b <path/to/sequence/database/file>, --seqdb <path/to/sequence/database/file>
						path of alignement file contrain by structure
						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form, fasta2dict

############# Fonction ####################
def filtreHit(files,dico_DB):
	"""
	"""
	nb = 0
	f = open(files,'r')
	lines = f.readlines()
	f.close()
	f = open('%s_filtred'%files,'w')
	files = open('%s_remove'%files,'w')
	goodAlignement = False
	for line in lines :
		if line[0]== '#':
			if line[0:2] == '#=' :
				if goodAlignement == True :
					f.write(line)
			else :
				f.write(line)
		elif line == '\n' :
			f.write(line)
		elif line[0:2] == '//' :
			f.write(line)

		else :
			name = line.split()[0]
			seq = dico_DB[name.split('/')[0]]
			index = 0
			goodAlignement = False
			for aa in seq :
				if aa == 'C' :
					zone = seq[index+33:index+49]
					if 'C' in zone :		
						f.write(line)
						goodAlignement = True
						break	
				index += 1
			
			if goodAlignement == False :
				nb +=1
				files.write(line)

	f.close()
	return nb

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to search a HMM profil for search EffecteurMAX. This program used HMMbuild et HMMsearch.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display HmmerPiepline version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-a', '--alignementFile',type = str, required=True, dest = 'file', help = 'Path of alignement file contrain by structure')
	filesreq.add_argument('-b', '--seqbd',type = str, required=True, dest = 'db', help = 'Path of sequence database')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')
	

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	alignement = os.path.abspath(args.file)
	db = os.path.abspath(args.db)
	outDir= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	outDir = verifDir(outDir)
	name_directory = [outDir,outDir]
	createDir(name_directory)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in HmmerPiepline (Version " + version + ")         ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')


############## main ####################################
	#os.system('module load bioinfo/hmmer/3.1b2')
	profilZero = '%sprofil_0'%outDir
	os.system('hmmbuild --amino %s %s > trash '%(profilZero,alignement))
	AlignementZero = '%salignement_0'%outDir
	os.system('hmmsearch -A %s --max --nonull2 %s %s > trash'%(AlignementZero,profilZero,db))
	dico_fasta = fasta2dict(db)
	nb = filtreHit(AlignementZero,dico_fasta)
	print('Alignement_0 : %s hits retirés'%(nb)) 
	newAlignement = AlignementZero
	i = 0
	while True:
		i += 1
		profil = '%sprofil_%s'%(outDir,str(i))
		os.system('hmmbuild --amino %s %s_filtred > trash'%(profil,newAlignement))
		f = open(newAlignement+'_filtred','r')
		linesOldAlignement = f.readlines()
		f.close()
		newAlignement = '%salignement_%s'%(outDir,str(i))
		os.system('hmmsearch -A %s --max --nonull2 -E 1e-3 %s %s > trash'%(newAlignement,profil,db))
		nb = filtreHit(newAlignement,dico_fasta)
		f = open(newAlignement+'_filtred','r')
		linesNewAlignement =  f.readlines()
		f.close()
		print('%s : %s hits retirés'%(newAlignement,nb)) 
		if linesOldAlignement == linesNewAlignement :
			print(i)
			break
		
		
	
			
	






















