#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package script_pangenome.py
# @author Florian Charriat

"""
	The script_pangenome script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/07/2018
	:version: 0.1

	Script description
	------------------

	This program is used to construct a pan génome. This script create a core genome File and a accessory genome file which contain only the group of orthologue of Effecteur MAX.
	
	Example
	-------

	>>> script_pangenome.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-a <path/to/alingment/HMM/file>, --align <path/to/fasta/file>
						path of the alignement HMM file of HmmerPipeline.py
		- \-db <path/to/blast/database>, --database <path/to/prosite.dat/file>
						path of the database of all genome which must be contain the database crate with makeblastbd.
		- \-f <path/to/secretome/fasta/file/>, --fasta <path/to/prosite.dat/file>
						path of the fasta file of all secretome
		- \-g                             --orthoGroups
						path of the directory which contain OrthoFinder result.
						
		- \-p                             --phenotype
						path of the phenotype file. Please, change the use the model for create your phenotype file.						
					
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
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
#	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to predict secretome with SignalP, TagetP, Phobius.
#	
#	This program uses the comparaisonSecretome script to retrieve and compare information from the secretome prediction tools.
#	
#	And the selection_TMHMM script to select, from TMHMM ouput, only protein with no TH domain or only one TH in 60 first aa.
#	
#	This program uses also the eliminateREmotif script to eliminate, from ps_scan ouput, the protein with a RE retention motif.
#	
#	Please make sure that this script is present in the directory.''')
#	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
#'display script_pangenome version number and exit')
#
#
#	filesreq = parser.add_argument_group('Input mandatory infos for running')
#	filesreq.add_argument('-a', '--align',type = str, required=True, dest = 'alignFile', help = 'Path of the alignement HMM file of #HmmerPipeline.py)
#	filesreq.add_argument('-db', '--database',type = str, required=True, dest = 'db', help = 'Path of the database use with #HmmerPipeline.py. The database must be contain the database crate with makeblastbd and the fasta file.')
#	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of the fasta file of all secretome')
#	filesreq.add_argument('-g', '--orthoGroups',type = str, required=True, dest = 'orthoGroupe', help = 'Path of the directory which #contain OrthoFinder result.')
#	filesreq.add_argument('-p', '--phenotype', action='store_true', dest = 'force', help = 'Path of the phenotype file. Please, use #the model for create your phenotype file')
#	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')


	
######### Recuperation arguments ###########
#	args = parser.parse_args()
#	pathAlignement = os.path.abspath(args.alignFile)
#	DB_blast= os.path.abspath(args.db)
#	pathSecretome = os.path.abspath(args.fasta)
#	prosite_dat =  os.path.abspath(args.prositePath)
#	verifFichier(prosite_dat)
	
	
########### Gestion directory ##############
#	verifFichier(fasta_file)
#	outDir = verifDir(outDir)
#	bash = outDir+'script_bash'
#	name_directory = [outDir,outDir+'error_files', outDir+'out_files',bash]
#	createDir(name_directory)


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("      Welcome in secretome_Pipeline (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')


	
	############### Params ############
	#pathAlignement = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v2/result_without_2_hits/alignement_7_filtred'
	#pathAlignement = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v3/2_result_isolat_one_by_one/merge_result'
	pathAlignement ='/homedir/gladieux/work/magMax_project/6_EffecteurMax/6_test_20_per_file/merge_first_alignement.fasta'
	#pathAlignement = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v1/test_profil_jerome'
	#pathAlignement = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v1/test_pipeline_on_Jerone_alignement/alignement_13_filtred'
	pathSecretome = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v1/secretome_DB.fasta'
	#pathOrtho = '/homedir/gladieux/work/magMax_project/4_Orthologie/0_rawdata/Results_Jul10_1/Orthogroups.csv'
	#pathOrtho = '/homedir/gladieux/work/magMax_project/4_Orthologie/1_selectRawdata/Results_msa__Jul25/Orthogroups.csv'
	pathOrtho = '/homedir/gladieux/work/magMax_project/4_Orthologie/4_test_MCL_2/0_rawdata/Results_MCL2_Jul30_1/Orthogroups_only_one_transcript_by_groups.txt'
	pathFastaQuery = '/homedir/gladieux/work/magMax_project/4_Orthologie/6_effecteurMax_pangenome/tmp/'
	DB_blast ='/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/0_DB/All_genome_isolat'
	output_blast = '/homedir/gladieux/work/magMax_project/4_Orthologie/6_effecteurMax_pangenome/tmp/'
	PathDataH = '/homedir/gladieux/work/magMax_project/4_Orthologie/6_effecteurMax_pangenome/hote_isolat.csv'
	#PathDataH = '/homedir/gladieux/work/magMax_project/4_Orthologie/6_effecteurMax_pangenome/Lineage_isolat.csv'
	PathOutput = '/homedir/gladieux/work/magMax_project/4_Orthologie/6_effecteurMax_pangenome/data-final/'
	pathScript = sys.path[0]+'/Report_pangenome.Rmd'
	#PathCount = '/homedir/gladieux/work/magMax_project/4_Orthologie/0_rawdata/Results_Jul10_1/Orthogroups.GeneCount.csv'
	#PathCount = '/homedir/gladieux/work/magMax_project/4_Orthologie/1_selectRawdata/Results_msa__Jul25/Orthogroups.GeneCount.csv'
	PathCount = '/homedir/gladieux/work/magMax_project/4_Orthologie/4_test_MCL_2/0_rawdata/Results_MCL2_Jul30_1/Orthogroups.GeneCount_only_one_transcript_by_groups.csv'


	
	############### Création file sequence Groupe + file Alignement (effecteur MAX) ###############################
	query = True # Si les fichiers sont déjà crée, query = True

	if query == True:
		print(form("1.  Recherche des groupe d'othologue d'effecteur MAX",'green','bold')+' (Déja réalisé)\n')
	if query == False : 
		print("1.  Recherche des groupe d'othologue d'effecteur MAX")
		############ Création liste Effecteurs MAX trouvé par recherche HMM ###############
		f = open(pathAlignement,'r')
		lines = f.readlines()
		f.close()
		effecteurMax = []
		for line in lines :
			if line[0] != '#' and line[0] != '\n' :
				effecteurMax.append(line.split('/')[0])
		############ Création dictionnaire de tous le secretome ############################
		dico_secretome = fasta2dict(pathSecretome)

		############### Création file sequence Groupe + file Alignement (effecteur MAX) ###############################
		f = open(pathOrtho,'r')
		lines = f.readlines()
		f.close()
		lines = lines[1:]
		for line in lines :
			groupe = line.split(':')[0]
			count = line.split()[1:]
			for elt in count :
				if elt in effecteurMax :
					fileName = pathFastaQuery+groupe+'.fasta'
					f = open(fileName,'w')
					sequence = Seq(str(dico_secretome[elt].seq) + '*')
					record = SeqRecord(sequence,id=str(groupe),name=str(groupe), description= 'length : '+ str(len(sequence)))
					SeqIO.write(record,f, "fasta")
					f.close()
					output = output_blast + groupe+'_blast.txt'
					os.system('tblastn -query %s -db %s -evalue 1e-4 > %s'%(fileName,DB_blast,output))	
					break
		f.close()	

	#################### Récupération info données blast ##############################################	
	print(form('2.  Correction des Annotations a partir des blasts effectués\n','green','bold'))
	listeGroupeEffecteur = []
	dico_groupe = {}
	for blast_file in os.listdir(output_blast):
		if blast_file.endswith('_blast.txt'):
			listeGroupeEffecteur.append(blast_file.replace('_blast.txt',''))
			startFile = True
			pathFile = output_blast+blast_file
			f = open(pathFile)
			lines = f.readlines()
			f.close()
			pathResult = output_blast+blast_file.replace('_blast.txt','_gene.txt')
			result = open(pathResult,'w')
			result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('name','Num_scaffold','start','end','sens','PcIdent','couverture'))
			liste = []
			for line in lines :
				if line[0:6] == 'Query=' :
					lenQuery = line.split()[-1]
				
				elif line[0] == '>' and startFile == True:
					#print('Start '+blast_file)
					lineSplit = line.split()
					name = lineSplit[1].split('_')[0]
					Num_scaffold = lineSplit[1].split('_')[-1]
					start = True
					startFile = False
				
				elif line[0] == '>' :
					if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
						start,end,sens = functionSens(pos1,pos2)
						result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name,Num_scaffold,start,end,sens,PcIdent,couverture))
						liste.append(name)
					lineSplit = line.split()
					name = lineSplit[1].split('_')[0]
					Num_scaffold = lineSplit[1].split('_')[-1]
					start = True
				elif 'Identities' in line :
					lineSplit = line.split()
					PcIdent = int(lineSplit[3].replace('(','').replace('%),',''))
					couverture = int(lineSplit[2].split('/')[-1])/int(lenQuery) *100
				elif 'Sbjct' in line and start == True:
					lineSplit = line.split()
					pos1 = lineSplit[1]
					startAA = lineSplit[2][0]
					start = False
				elif 'Sbjct' in line and start == False:
					lineSplit = line.split()
					pos2 = lineSplit[1]
					endAA = lineSplit[2][-1]
			if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
				start,end,sens = functionSens(pos1,pos2)
				result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name,Num_scaffold,start,end,sens,PcIdent,couverture))
				liste.append(name)
			result.close()
			dico_groupe[blast_file.replace('_blast.txt','')] = liste	
				
			
	############ Création dictionnaire des counts des groupes othologue ################
	print(form("3.  Création du Core genome et des groupes orthologues d'Effecteur MAX\n",'green','bold'))
	f = open(PathCount,'r')
	lines = f.readlines()
	f.close()
	line0 = 'Groupe\t'+ '\t'.join(lines[0].replace('_protein','').replace('Total','').split())+'\n'
	core = open(PathOutput+'coreGenome.csv','w')
	core.write(line0)
	max = open(PathOutput+'MaxGenome.csv','w')
	max.write(line0)
	entete = lines[0].split()
	lines = lines[1:]
	dico_count = {}
	for line in lines :
		groupe = line.split()[0]
		count = line.split()[1:-1]
		liste_other = []
		if '0' not in count :
			core.write(line)
		if groupe in listeGroupeEffecteur :
			listeIndex = indexEgale(count,'0')
			for i in listeIndex :
				if entete[i].replace('_protein','') in dico_groupe[groupe] :
					count[i] = str(dico_groupe[groupe].count(entete[i].replace('_protein','')))+'*'
			max.write('%s\t%s\n'%(groupe,'\t'.join(count)))
		dico_count[groupe]=count
	core.close()	
	max.close()	
	print(form('4.  Création du  rapport\n','green','bold'))
	os.system("Rscript -e 'rmarkdown::render(%s, output_file=%s, quiet=TRUE, params = list(data = %s, dataH = %s))'"%('"'+pathScript+'"','"'+PathOutput+'Report.html'+'"','"'+PathOutput+'MaxGenome.csv"','"'+PathDataH+'"'))
	print("Rscript -e 'rmarkdown::render(%s, output_file=%s, quiet=TRUE, params = list(data = %s, dataH = %s))'"%('"'+pathScript+'"','"'+PathOutput+'Report.html'+'"','"'+PathOutput+'MaxGenome.csv"','"'+PathDataH+'"'))


############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))







