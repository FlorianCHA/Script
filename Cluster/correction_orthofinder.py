#!/usr/local/bioinfo/python/3.4.3_build2/bin/pythonqs
# -*- coding: utf-8 -*-
# @package comparaisonSecretome.py
# @author Florian Charriat


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens

##########Path ###################

pathFasta = '/homedir/gladieux/work/magMax_project/6_EffecteurMax/3_psiblastAllprot/DB/all_protein.fasta'
pathDBGenome = ''
pathOrthologue = '/homedir/gladieux/work/magMax_project/4_Orthologie/1_selectRawdata/Results_msa__Jul25/'
ouput = '/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/'
tmp = '/homedir/gladieux/work/magMax_project/4_Orthologie/3_correction/tmp/'




########################## Création dico absence/présence ####################################################################
f = open(pathOrthologue+'Orthogroups.GeneCount.csv','r')
lines = f.readlines()
f.close() 
entete = lines[0]
listeProtein = []
for elt in entete :
	if 'Total' not in elt :
		listeProtein.append(elt.replace('_protein','')

OG = lines[1:len(lines)]
dico_OG = {}
for line in lines :
	OG = line.split('\t')[0]
	count = line.split('\t')[1:-1]
	if '0' not in count :
		dico_OG[OG] = listeProtein
	else :
		listeIndex = indexDif(count,'0')
		listeCount = []
		for i in listeIndex :
			listeCount.append(listeProtein[i])
		dico_OG[OG] = listeCount
			

############### Correction ############################################################################################

f = open(pathOrthologue+'Orthogroups.txt','r')
lines = f.readlines()
f.close()
lines = lines[1:len(lines)]
dico_secretome = fasta2dict(pathFasta)
dico_groupe = {}
for line in lines :
	listeLine = line.split()
	OG = listeLine[0].replace(':','')
	Idseq = listeLine[1].replace(' ','')
	fasta_output = tmp+ OG + '.fasta'
	if len(dico_OG[OG]) != len(listeProtein):
		dico_groupe[OG] = []
		f = open(fasta_output,'w')
		sequence = Seq(str(dico_secretome[Idseq].seq) + '*')
		record = SeqRecord(sequence,id=str(OG),name=str(OG), description= 'length : '+ str(len(sequence)))
		SeqIO.write(record,f, "fasta")
		f.close()
		blast_output = tmp + OG +'_blast.txt'
		print('tblastn -query %s -db %s -evalue 1e-4 > %s'%(fasta_output,pathDBGenome,blast_output))	
		f = open(blast_output)
		lines1 = f.readlines()
		f.close()
		pathResult = output+OG+'_gene.txt'
		result = open(pathResult,'w')
		result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('name','Num_scaffold','start','end','sens','PcIdent','couverture'))
		liste = []
		for line in lines1 :
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
			if name not in dico_OG[OG] :
				start,end,sens = functionSens(pos1,pos2)
				result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name,Num_scaffold,start,end,sens,PcIdent,couverture))
				dico_groupe[OG].append('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end))
		result.close()
	
f = open(pathOrthologue+'Orthogroups.csv','r')
lines = f.readlines()
f.close()
result = open(pathOrthologue+'Orthogroups_correction.csv','w')
for line in lines :
	OG = line.split('\t')[0]
	newline = line.replace('\n',',') + ', '.join(dico_groupe[OG])
	result.write(newline+'\n')
	
result.close()

f = open(pathOrthologue+'Orthogroups.txt','r')
lines = f.readlines()
f.close()
result = open(pathOrthologue+'Orthogroups_correction.txt','w')
for line in lines :
	OG = line.split(':')[0]
	newline = line.replace('\n',' ') + ' '.join(dico_groupe[OG])
	result.write(newline+'\n')
result.close()

f = open(pathOrthologue+'Orthogroups.GeneCount.csv','r')
lines = f.readlines()
f.close()
result = open(pathOrthologue+'Orthogroups.GeneCount_correction.txt','w')
for line in lines :
	OG = line.split('\t')[0]
	lineSplit = line.split('\t')
	for elt in dico_groupe[OG]:
		listeI = indexDif(listeProtein,elt.split('_')[0])
		for i in listeI :
			lineSplit[i] = '1'
	
	newline = '\t'.join(lineSplit)
	result.write(newline+'\n')
result.close()
	
	
	
	
	
	
	
	
