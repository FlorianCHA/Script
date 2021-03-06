#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package OG_phylo.py
# @author Florian Charriat

"""
	The ABYSS_launch script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 21/08/2018
	:version: 0.1

	Script description
	------------------

	This script was created for a special case (phylogenetic OG). So they have no documentation. Please contact the author for more information
"""


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from progress.bar import ChargingBar

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens


################## PATH ##################################################
output_tmp = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/tmp/'
MGG_fasta = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/1_MGG/Magnaporthe_oryzae.MG8.pep.all.fasta'
pathDBGenome = '/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/0_DB/All_genome_isolat'
pathOrthologue = '/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/'

###############################################################################################################
################## Récupération des groupes nescessaire à la phylogénétique ###################################
###############################################################################################################
f = open(pathOrthologue+'76Mo_2.txt','r')
lines = f.readlines()
f.close()
OG_phylo = []
print('')
bar = ChargingBar('Récupération des %s groupes Orthologue'%len(lines), max=len(lines), suffix='%(percent)d%%')
for OG in lines :
	OG = OG.replace('_Orthologue_macse_NT.fasta\n','')
	OG_phylo.append(OG)
	bar.next()
bar.finish()	

###############################################################################################################
#################### récupération des counts des OG ###########################################################
###############################################################################################################

f = open(pathOrthologue+'Orthogroups_1.GeneCount_new.csv','r')
lines = f.readlines()
f.close()
entete = lines[0].split()
lines = lines[1:len(lines)]
listeProtein = []
for elt in entete :
	if 'Total' not in elt :
		listeProtein.append(elt.replace('_protein',''))
dico_count = {}
dico_OG = {}
print('')

###############################################################################################################
##################### Récupération des OGs séléctionnés qui ne sont pas dans le core genome ###################
###############################################################################################################

bar = ChargingBar('Lecture du fichier GeneCount.csv', max=len(lines), suffix='%(percent)d%%')
for line in lines :
	listOG = line.split()
	OG = listOG[0]
	dico_count[OG] = listOG[1:len(listOG)]
	count = line.split('\t')[1:-1]
	if '0' not in count :
		dico_OG[OG] = listeProtein
	else :
		listeIndex = indexDif(count,'0')
		listeCount = []
		for i in listeIndex :
			listeCount.append(listeProtein[i])
		dico_OG[OG] = listeCount
	bar.next()
bar.finish()	


######################################################################################################################
####################### Création des fichiers fasta pour chaque OG ###################################################
######################################################################################################################

f = open('/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/Orthogroups_1.txt','r')
lines = f.readlines()
f.close()


dico_MGG = fasta2dict(MGG_fasta)
nb = 0
dico_OGvsMGG = {} 
OG_76MO = []
listeOG = []
result = open(pathOrthologue+'Orthogroups_phylo.txt','w')
print('')
bar = ChargingBar('Création des fichiers fasta de tous les OG séléctionnées', max=len(OG_phylo), suffix='%(percent)d%%')
for line in lines :
	for OG in OG_phylo :
		if OG in line :
			dico_OGvsMGG[OG] = line.split(':')[0]
			listeOG.append(OG)
			newLine = line.replace(OG,'')
			newline = newLine.replace(line.split(':')[0],OG)
			result.write(newline)
			
			if line.split(':')[0] not in dico_count.keys() or '0' in dico_count[line.split(':')[0]]:
				OG_output = output_tmp +'fasta/'+OG+'.fasta'
				OG_76MO.append(OG)
				f = open(OG_output,'w')
				sequence = Seq(str(dico_MGG[OG].seq) + '*')
				record = SeqRecord(sequence,id=str(OG),name=str(OG), description= 'length : '+ str(len(sequence)))
				SeqIO.write(record,f, "fasta")
			bar.next()
			nb += 1
			break

	if nb == len(OG_phylo) :
		break
	
bar.finish()
result.close()

############################################################################################################################
################################### Création des counts et traitement de ces dernier pour la correction ####################
############################################################################################################################

f = open(pathOrthologue+'/Orthogroups_1.GeneCount_new.csv','r')
lines = f.readlines()
f.close()
entete = lines[0].split()
print('')
bar = ChargingBar('Création du fichier et du dico des counts des OG séléctionnés', max=len(listeOG), suffix='%(percent)d%%')
lines = lines[1:len(lines)]
result = open(pathOrthologue+'Orthogroups_phylo.GeneCount.txt','w')
dico_OG = {}
dico_count = {}
for line in lines :
	for OGp in listeOG :
		if dico_OGvsMGG[OGp] in line :
			listOG = line.split()
			newline = line.replace(dico_OGvsMGG[OGp],OGp)
			result.write(newline)
			OG = listOG[0]
			count = line.split('\t')[1:-1]
			dico_count[OGp] = listOG[1:len(listOG)]
			if '0' not in count :
				dico_OG[OGp] = listeProtein
			else :
				listeIndex = indexDif(count,'0')
				listeCount = []
				for i in listeIndex :
					listeCount.append(listeProtein[i])
				dico_OG[OGp] = listeCount
			bar.next()
			break

bar.finish()		
result.close()
###########################################################################################################################
############################## Partie lancement Blast pour correction #####################################################
###########################################################################################################################

blast = False 
allscript = open(pathOrthologue+'tmp/LaunchBlast.sh','w')
if blast :
	for OG in OG_76MO :
		fasta_output = output_tmp +'fasta/'+OG+'.fasta'
		OG_output = output_tmp +'script/'+OG+'.sh'
		blast_output = output_tmp +'blast/'+OG+'.txt'
		outputfile =  output_tmp +'output_file/'+OG
		script = open(OG_output,'w')
		script.write('#$ -o %s.o\n#$ -e %s.e\n'%(outputfile,outputfile))
		script.write('module load bioinfo/ncbi-blast/2.6.0\n')
		script.write('tblastn -query %s -db %s -evalue 1e-4 > %s'%(fasta_output,pathDBGenome,blast_output))	
		script.close()
		os.system('chmod 755 %sscript/*'%output_tmp)
		allscript.write('qsub -q normal.q -V -cwd %s\n'%OG_output)
	
else :
	print('\nBlast already done')

allscript.close()

#################################################################################################################
######################################### Parsage des résultat blast + ##########################################
#################################################      +       ##################################################
######################################### Création fasta des correction #########################################
#################################################################################################################

result = open(output_tmp+'gene_add.txt','w')
result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('name','Num_scaffold','start','end','sens','PcIdent','couverture'))
seq = ''
dico_groupe = {}
dico_fasta = {}
nbC = 0
print('')
bar = ChargingBar('Correction des OG séléctionnés ', max=len(OG_76MO), suffix='%(percent)d%%')
for OG in OG_76MO :
	dico_groupe[OG] = []
	dico_fasta[OG] = []
	blast_output = output_tmp +'blast/'+OG+'.txt'
	f = open(blast_output,'r')
	lines = f.readlines()
	f.close()
	f = open('%scorrection/%s.fasta'%(output_tmp,OG),'w')
	liste = []
	startFile = True
	for line in lines :
		if line[0:6] == 'Query=' :
			lenQuery = line.split()[-1]

		elif line[0] == '>' and startFile == True:
			lineSplit = line.split()
			name = lineSplit[1].split('_')[0]
			Num_scaffold = lineSplit[1].split('_')[-1]
			start = True
			startFile = False
		elif line[0] == '>' :
			if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
				if name not in dico_OG[OG] :
					start,end,sens = functionSens(pos1,pos2)
					result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name,Num_scaffold,start,end,sens,PcIdent,couverture))
					dico_groupe[OG].append('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end))
					dico_fasta[OG].append(['%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end),seq])
					sequence = Seq(seq)
					record = SeqRecord(sequence,id=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)),name=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)), description= 'length : '+ str(len(sequence)))
					SeqIO.write(record,f, "fasta")
					nbC += 1

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
			seq = lineSplit[2]
			start = False
		elif 'Sbjct' in line and start == False:
			lineSplit = line.split()
			pos2 = lineSplit[1]
			endAA = lineSplit[2][-1]
			seq = seq +  lineSplit[2]
		
	if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
		if name not in dico_OG[OG] :
			start,end,sens = functionSens(pos1,pos2)
			result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name,Num_scaffold,start,end,sens,PcIdent,couverture))
			dico_groupe[OG].append('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end))
	f.close()
	bar.next()
bar.finish()
result.close()

print('\nNombre de correction effectué sur les %s Groupes Orthologue  : %s\n'%(len(OG_76MO),nbC))

##################################################################################################################
##################################### Ecriture nouveau fichier avec correction ###################################
##################################################################################################################

f = open(pathOrthologue +'Orthogroups_phylo.txt','r')
lines = f.readlines()
f.close()
lines = lines[1:len(lines)]
result = open(output_tmp+'Orthogroups_correction.txt','w')
bar = ChargingBar('Création Orthogroups_correction.txt', max=len(lines), suffix='%(percent)d%%')
for line in lines :
	OG = line.split(':')[0]
	if OG in OG_76MO and len(dico_groupe[OG]) != 0 :
		newline = line.replace('\n',' ') + ' '.join(dico_groupe[OG])
	else :
		newline = line.replace('\n',' ')
	result.write(newline+'\n')
	bar.next()
bar.finish()
result.close()

f = open(pathOrthologue+'Orthogroups_phylo.GeneCount.txt','r')
lines = f.readlines()
f.close()
lines = lines[1:len(lines)]
result = open(output_tmp +'Orthogroups.GeneCount_correction.csv','w')
bar = ChargingBar('Création Orthogroups.GeneCount_correction.csv', max=len(lines), suffix='%(percent)d%%')
for line in lines :
	OG = line.split('\t')[0]
	lineSplit = line.split('\t')[1:len(line.split('\t'))]
	if OG in OG_76MO and len(dico_groupe[OG]) != 0 :
		for elt in dico_groupe[OG]:
			listeI = indexEgale(listeProtein,elt.split('_')[0])
			for i in listeI :
				lineSplit[i] = '1*'
	
		newline = '\t'.join(lineSplit)
	else :
		newline = line
	result.write(newline+'\n')
	bar.next()
bar.finish()
result.close()

