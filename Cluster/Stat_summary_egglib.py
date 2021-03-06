#!/usr/local/bioinfo/python/2.7.9/bin/python
# -*- coding: utf-8 -*-
# @package summary_stat_egglib.py
# @author Florian CHARRIAT

"""
	The summary_stat_egglib script
	===========================
	:author: Florian CHARRIAT & Pierre Gladieux
	:contact: florian.charriat@inra.fr
	:date: 08/07/2018
	:version: 0.1
	Script description
	------------------
	This Programme is used to correct vcf file with multi sample create with show-snp
	Example
	-------
	>>> summary_stat_egglib.py -f multiSample.vcf -o output.vcf -d show-snp.result


	Help Programm
	-------------
	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display summary_stat_egglib.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/show-snp/result>, --file <path/to/show-snp/result>
						path of the vcf file with multi sample
		- \-d <path/to/show-snp/directory>, --directory <path/to/show-snp/directory>
						path of the directory which contains all show-snp file without -I and -C option.
		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file ( format vcf)

"""


##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,is_number

## Python modules
import argparse, egglib
from time import gmtime, strftime
# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to calculate diversity statistics between groups with the Egglib package''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'directory', help = 'path of the directory which contains all the fasta file to proceed')
	filesreq.add_argument('-g', '--groups',type = str, required=True, dest = 'groups', help = 'path of the file which contain per line : Name of Isolat\tGroups')
	filesreq.add_argument('--focal_groups', type=str,nargs='+', default = 'None', dest='target',
						 help='One or more groups to analyse')
	filesreq.add_argument('--min_length', type=str, required=True, dest='minLength',
						 help='the minimum length of sequence to analyse')
	filesreq.add_argument('--min_size', type=str, required=True, dest='minSize',
						 help='the minimum size of exploitable samples')
	filesreq.add_argument('--translatorX', action = 'store_true', dest='translatorX',
						 help="If in your fasta file, your sequence havn't the same length please use this option for use transltorX option for align all sequence")
	filesreq.add_argument('--MISS', type=str, required=True, dest='miss',
						 help='maximum proportion of missing data (if there are more missing data at a site, the site is ignored)')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of the output file ( format vcf)')


######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.directory)
	groups = os.path.abspath(args.groups)
	target = args.target
	translatorX = args.translatorX
	MISS = args.miss
	minLength = args.minLength
	minSize = args.minSize
	output = os.path.abspath(args.output)
########### Gestion directory ##############
	directory = verifDir(directory,True)
	output = verifDir(output)
	createDir([output,output+'/translatorX',output+'/fastaByGroup'])
############### start message ########################

	print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
	print(
			"\t" + form("|", 'yellow', 'bold') + form("       Welcome in summary_stat_egglib (Version " + version + ")     ",
													  type='bold') + form("|", 'yellow', 'bold'))
	print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
	entete_file = 'Name\tNb samples\tNb analysed sites\tAverage nb ofsamples used\tNb segregating sites\tNb haplotypes\tExpected heterozygosity\tWatterson’s estimator of theta\tNucleotide diversity\tTajima’s D\tNb Non-synonimous\tNb synonimous\tNSseg\tSseg\tNucleotide diversity of non-synonimous site (Pi)\tNucleotide diversity of synonimous site (Pi)\n'
	summary_stat = {}
	dico_group = {}
	OUTs = {}
	all_stats = ['n', 'lseff', 'nseff', 'S', 'K', 'He', 'thetaW', 'Pi', 'D', 'NSsites', 'Ssites', 'NSseg', 'Sseg',
				 'PiNS', 'PiS']

	print('Lecture du fichier : %s'% groups.split('/')[-1])
	with open(groups,'r') as group_file :
		for line in group_file :
			isolat,group = line.strip().split()
			if group in target or target == 'None':
				if group not in dico_group.keys() :
					# Initialisation du dico summary pour donner à la fin les moyennes par groupes
					summary_stat[group] = {}
					for stat in all_stats:
						summary_stat[group][stat] = []
					# Création dico contenant tous les isolats (value) par groupe (keys)
					dico_group[group] = [isolat]
					# Création d'un dico de output
					OUTs[group] = open(output + group + '.csv', 'w')
					OUTs[group].write(entete_file)
				else :
					dico_group[group].append(isolat)

	total = len([file for file in os.listdir(directory) if file.endswith('.fasta')])
	count = 0




	for file in os.listdir(directory):
		count +=1
		print('Traitement du fichier : %s (%s/%s)'%(file,count,total))

		#Ouverture du fichier fasta pour le mettre dans un dixo et le split par groupe
		try :
			dico_fasta_all = fasta2dict(directory + file)
		except :
			print('\t-The %s file have two sequence with the same \n'%file)
			continue
		file_exist = False
		i = 0
		haplotype = {}
		for gp in sorted(dico_group.keys()) :
			fasta_group = output+'fastaByGroup/' + '%s_%s.fasta'%(file.replace(".fasta",""),gp)
			with open(fasta_group,'w') as fasta_file :
				for elt in dico_fasta_all.keys() :
					name = elt.replace('Mo_','').split('_')[0]
					if name in dico_group[gp] :
						record = SeqRecord(dico_fasta_all[elt].seq, id=str(elt), name=str(elt), description='')
						SeqIO.write(record, fasta_file, "fasta")
						file_exist = True
			if file_exist == False :
				os.system('rm %s'%fasta_group)
				continue
			file_exist = False


			if translatorX :
				file_aln = output + 'translatorX/' + '%s_%s.fasta'%(file.replace(".fasta",""),gp)
				os.system('translatorx_vLocal.pl -i %s -o %s  > stdout.txt 2> stdout.txt' % (fasta_group,file_aln ))
				file_aln = file_aln + '.nt_ali.fasta'

			else :
				file_aln = directory + file

			dico_fasta = fasta2dict(file_aln)
			i += 1
			haplotype[gp] = [[elt,str(dico_fasta[elt].seq),i] for elt in dico_fasta.keys() if elt.replace('Mo_','').split('_')[0] in dico_group[gp]]
			liste = [len(dico_fasta[elt].seq) for elt in dico_fasta.keys() if elt.replace('Mo_','').split('_')[0] in dico_group[gp]]
			aln = egglib.Align.create(haplotype[gp])
			cs = egglib.stats.ComputeStats()
			cs.add_stats('lseff', 'nseff', 'S', 'K', 'He', 'thetaW', 'Pi', 'D')
			stat_order = ['lseff', 'nseff', 'S', 'K', 'He', 'thetaW', 'Pi', 'D']
			stats = cs.process_align(aln, max_missing=float(MISS))
			if stats['lseff']>=int(minLength) and stats['nseff']>=int(minSize):
				if stats['Pi'] != None:
					Pi = stats['Pi'] / float(stats['lseff'])
				else:
					Pi = 0
				if stats['thetaW'] != 'None':
					thetaW = stats['thetaW'] / float(stats['lseff'])
				else:
					thetaW = 0
				stats_treated = stats
				# for stat in stats:
				# 	if stat == 'thetaW':
				# 		stats_treated['thetaW'] = thetaW
				# 	elif stat == 'Pi':
				# 		stats_treated['Pi'] = Pi
				# 	else:
				# 		stats_treated[stat] = stats[stat]
				rf = egglib.tools.ReadingFrame([(0, stats['lseff'] - 1)])
				cd = egglib.stats.CodingDiversity(aln, frame=rf) # Probleme Rf prend le 0 de l'alignement donc prend pas toute le séquence mais les GAP aussi (donc manque des séquence + ajout plein de gap)
				align_S = cd.mk_align_S()
				align_NS = cd.mk_align_NS()
				numS = cd.num_sites_S
				numNS = cd.num_sites_NS
				codon_filter = egglib.stats.Filter(rng=(0, 63), missing=64)
				cs = egglib.stats.ComputeStats()
				cs.add_stats('Pi', 'lseff')
				statsS = cs.process_align(align_S, filtr=codon_filter, max_missing=float(MISS))
				statsNS = cs.process_align(align_NS, filtr=codon_filter, max_missing=float(MISS))
				if statsS['Pi'] != None:
					PiS = statsS['Pi'] / numS
				else:
					PiS = 0
				if statsNS['Pi'] != None:
					PiNS = statsNS['Pi'] / numNS
				else:
					PiNS = 0
				OUTs[gp].write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(file.replace('.fasta', ''),aln.ns,'\t'.join([str(stats_treated[x]) for x in stat_order]),round(cd.num_sites_NS, 2),round(cd.num_sites_S, 2),cd.num_pol_NS,cd.num_pol_S,PiNS,PiS))


				# Création des données pour donner la moyenne des statistiques par groupe
				summary_stat[gp]['n'].append(aln.ns)
				for stat in stat_order :
					if stats_treated[stat] != None:
						summary_stat[gp][stat].append(stats_treated[stat])
					else:
						summary_stat[gp][stat].append(0)
				summary_stat[gp]['NSsites'].append(round(cd.num_sites_NS, 2))
				summary_stat[gp]['Ssites'].append(round(cd.num_sites_S, 2))
				summary_stat[gp]['NSseg'].append(cd.num_pol_NS)
				summary_stat[gp]['Sseg'].append(cd.num_pol_S)
				summary_stat[gp]['PiNS'].append(PiNS)
				summary_stat[gp]['PiS'].append(PiS)

	entete_file = 'Name\tNb sequences used\tNb samples\tNb analysed sites\tAverage nb ofsamples used\tNb segregating sites\tNb haplotypes\tExpected heterozygosity\tWatterson’s estimator of theta\tNucleotide diversity\tTajima’s D\tNb Non-synonimous\tNb synonimous\tNSseg\tSseg\tNucleotide diversity of non-synonimous site (Pi)\tNucleotide diversity of synonimous site (Pi)\n'
	with open(output+'summary_stat.csv','w') as summary_file :
		summary_file.write(entete_file)
		for gp in OUTs.keys():
			OUTs[gp].close()
			stats_mean = [str(len(summary_stat[gp]['n']))]
			for stat in all_stats :
				liste_stat = [float(x) for x in summary_stat[gp][stat]]
				if len(liste_stat) != 0 :
					stats_mean.append(str(sum(liste_stat) / len(liste_stat)))
				else :
					stats_mean.append('No data')
			summary_file.write('%s\t%s\n'%(gp,'\t'.join(stats_mean)))




	############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))