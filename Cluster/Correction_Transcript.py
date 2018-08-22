#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package Correction_Transcript.py
# @author Florian Charriat

########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human	
from progress.bar import ChargingBar

f = open('/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/Orthogroups_1.txt','r')
lines = f.readlines()
f.close()


result = open('/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/Orthogroups_1_new.txt','w')
dico = {}
dico_count = {}

print('')

bar = ChargingBar('Processing for Orthogroups.txt ', max=len(lines), suffix='%(percent)d%%')

for line in lines :
	protein_done = []
	listOG = line.split()
	OG = listOG[0]
	dico[OG] = []
	dico_count[OG.replace(':','')] = []
	listProtein = listOG[1:len(listOG)]
	for protein in listProtein :
		if 'T0' in protein :
			protein_done.append(protein[0:-2])
			dico[OG].append(protein)
			if 'MGG' in protein :
				dico_count[OG.replace(':','')].append('MGG')
				
			else : 
				dico_count[OG.replace(':','')].append(protein.split('_')[1])

	result.write('%s %s\n'%(OG, ' '.join(dico[OG])))
	bar.next()
bar.finish()	
result.close()


f = open('/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/Orthogroups_1.GeneCount.csv','r')
lines = f.readlines()
f.close()
print('')
entete = []
for elt in lines[0].split('\t') :
	if 'Magnaporthe_oryzae.MG8' in elt :
		entete.append('MGG')
	else :
		entete.append(elt.split('_')[0])

result = open('/homedir/gladieux/work/magMax_project/4_Orthologie/5_test_without_transcript/0_rawdata/Results_MCL2_T0_Aug06/Orthogroups_1.GeneCount_new.csv','w')
result.write('\t'.join(entete))
print('')
bar = ChargingBar('Processing for Orthogroups.GeneCount.csv ', max=len(dico_count.keys()), suffix='%(percent)d%%')
for OG in sorted(dico_count.keys(), key=sort_human):
	listeCount = []
	Total = 0
	for elt in entete[1:(len(entete)-1)] :
		if elt == 'Magnaporthe-Wheatblast-12-1-205' :
			listeCount.append(str(dico_count[OG].count('Magnaporphe-Wheatblast-12-1-205')))
		else :
			listeCount.append(str(dico_count[OG].count(elt)))
	result.write('%s\t%s\t%s\n'%(OG,'\t'.join(listeCount),len(dico_count[OG])))
	bar.next()
bar.finish()	
result.close()
	
	
	
	
	
