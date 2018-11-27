#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package makeEffectome.py
# @author Florian Charriat

"""
	The makeEffectome script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 15/10/2018
	:version: 0.1

	Script description
	------------------

	This program is used to create a fasta file which contain all nucléotidique sequence of effecteur MAX with 100bp before and after the sequence.

	Example
	-------

	>>> makeEffectome.py -f /homedir/user/work/fasta.file -p /homedir/user/work/getORF.result -o /homedir/user/work/result.gff

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display makeEffectome.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
						path of fasta that contains all the genome assembly of Isolat
						
		- \-p <path/to/fasta/file>, --prot <path/to/fasta/file>
						path of fasta that contains all the protein of Isolat.

		- \-hmm <path/to/hmmsearch/file>, --hmmsearch <path/to/fasta/file>
						path of hmmsearch result (stdout of effecteur MAX search)

		- \-o <path/to/output/file>, --output <path/to/output/file>
						path of the output file

"""


########## Module ###############
## Python modules
qimport argparse, os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from progress.bar import ChargingBar

#Import module_Flo
from module_Flo import verifDir, createDir , form,verifFichier,fasta2dict, openfile,is_number



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to create a fasta file which contain all nucléotidique sequence of effecteur MAX with 100bp before and after the sequence.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display makeEffectome version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of fasta file that contains all the genome assembly of Isolat')
	filesreq.add_argument('-p', '--prot',type = str, required=True, dest = 'protein', help = 'Path of fasta that contains all the protein of Isolat.')
	filesreq.add_argument('-hmm', '--hmmsearch',type = str, required=True, dest = 'hmm', help = 'Path of hmmsearch result (stdout of effecteur MAX search)')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	genome = os.path.abspath(args.fasta)
	prot = os.path.abspath(args.protein)
	hmm = os.path.abspath(args.hmm)
	output = os.path.abspath( args.output)

########### Gestion directory ##############
	verifDir(genome,True)
	verifFichier(prot)
	verifFichier(hmm)

################## Main ######################
	print('Récupération des séquences et informations de toute les protéines')
	lines = openfile(prot)
	dico_prot = {}
	print('')
	bar = ChargingBar('Recuperation des informations de chaque gène du fichier protéique : ', max=len(lines), suffix='%(percent)d%%')
	for line in lines :
		bar.next()
		if '>' in line :
			name = line.split('|')[0].replace('>','').strip()
			Isolat = name.split('_')[1]
			num_scaffold = line.split('_')[-3].strip()
			scaffold = 'Scaffold_%s'%(num_scaffold)
			position = line.split('_')[-2].split('|')[0].split(':')
			if is_number(position[0]) and len(position) == 2:
				start = int(position[0].strip())
				end = int(position[1].strip())
				length = line.split('length=')[-1].replace('\n','').strip()
				dico_prot[name] = [scaffold,start,end,length]

	bar.finish()	

	hmm_result = openfile(hmm)
	start_seq = False
	effecteur = []
	
	print('')
	print('Recuperation des sequences de chaque effecteur')
	for line in hmm_result :
		if '------' in line :
			start_seq = True

		elif start_seq == True and len(line.split()) < 9 :
			break
		
		elif start_seq == True:
			name = line.split()[8]
			if 'Magnaporphe-Wheatblast-12-1-205' not in name and 'MGG' not in name :
				effecteur.append(name)

	hmm_result = ''
	old_name = ''
	f = open(output,'w')
	for elt in sorted(effecteur):
		name = elt.split('_')[1]
		if name != old_name :
			print(name)
			fasta_genome = fasta2dict(genome+'/'+name+'.fasta')
			
		scaffold, start, end, length = dico_prot[elt]
		if int(length) < 300 :
			sequence = Seq(str(fasta_genome[scaffold].seq)[start-90:end+90])
			record = SeqRecord(sequence,id=str(elt),name=str(elt), description= 'length : '+ str(len(sequence)))
			SeqIO.write(record,f, "fasta")
		old_name = name

	f.close()



