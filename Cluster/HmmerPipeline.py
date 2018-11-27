#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package HmmerPiepline.py
# @author Florian Charriat

"""
	The HmmerPiepline script
	========================
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

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human

############# Fonction ####################
def filtreHit(files,dico_DB):
	"""
	"""
	dico = {}
	nb = 0
	nbU = 0
	nbG = 0
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
			aln = line.split()[1]
			index = 0
			goodAlignement = False
			for aa in seq :
				if aa == 'C' :
					zone = seq[index+33:index+50]
					if 'C' in zone :
						f.write(line)
						if name not in dico.keys():
							nbG += 1
							dico[name] = aln
						elif name in dico.keys():
							dico[name] = dico[name].replace('\n','') + aln
						goodAlignement = True
						break
				index += 1
		
			if goodAlignement == False :
				nb +=1
				files.write(line)


#	listeAln = []
#	for elt in sorted(dico.keys(), key=sort_human):
#		nbG += 1
#		if dico[elt] not in listeAln :
#			nbU += 1
#			listeAln.append(dico[elt])
#			f.write(elt+'\t'+dico[elt]+'\n')
#	f.write('//')
	f.close()
	files.close()
	return nbG,nb#,nbU

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
	name_directory = [outDir]
	createDir(name_directory)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in HmmerPiepline (Version " + version + ")         ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')


############## main ####################################
	os.system('module load bioinfo/hmmer/3.1b2')
	profilZero = '%sprofil_0'%outDir
	os.system('hmmbuild --amino %s %s > trash '%(profilZero,alignement))
	AlignementZero = '%salignement_0'%outDir
	newAlignement = AlignementZero
	print(form('Alignement n° 0\n','green','bold'))
	print(form('\t - Création du pattern 0 à partir de %s_filtred'% (newAlignement.split('/')[-1]),'white','bold'))	
	os.system('hmmsearch -A %s -E 1e-3 --max --nonull2 %s %s > trash'%(AlignementZero,profilZero,db))
	dico_fasta = fasta2dict(db)
	nbG,nb = filtreHit(AlignementZero,dico_fasta)
	print(form('\t - Alignement n° 0 filtré','white','bold'))	

	print(form('\t - %s : %s hits récupérés, %s hits retirés\n'%(newAlignement.split('/')[-1],nbG,nb),'white','bold'))	

	i = 0
	while True:
		i += 1
		print(form('Alignement n° %s\n'%i,'green','bold'))
		print(form('\t - Création du pattern %s a partir des %s hits du fichier %s_filtred'% (i,nbG,newAlignement.split('/')[-1]),'white','bold'))	
		profil = '%sprofil_%s'%(outDir,str(i))
		os.system('hmmbuild --amino %s %s_filtred > trash'%(profil,newAlignement))
		f = open(newAlignement+'_filtred','r')
		linesOldAlignement = f.readlines()
		f.close()
		newAlignement = '%salignement_%s'%(outDir,str(i))
		print(form('\t - Utilisation du profil n° %s sur le fichier %s'%(i,db.split('/')[-1]),'white','bold'))	
		os.system('hmmsearch -A %s --max --nonull2 -E 1e-3 %s %s > %s'%(newAlignement,profil,db,'%sresult_%s'%(outDir,str(i))))
		nbG,nb = filtreHit(newAlignement,dico_fasta)
		print(form('\t - Alignement n° %s filtré'%i,'white','bold'))	
		f = open(newAlignement+'_filtred','r')
		linesNewAlignement =  f.readlines()
		f.close()
		print(form('\t - %s : %s hits récupérés, %s hits retirés\n'%(newAlignement.split('/')[-1],nbG,nb),'white','bold'))	
		if linesOldAlignement == linesNewAlignement or len(linesOldAlignement) > len(linesNewAlignement):
			print('Il y a eu %s itération'%i)
			break
		
		
	
			
	


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




















