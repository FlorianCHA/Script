#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package QualityAssemblage.py
# @author Florian Charriat

"""
	The QualityAssemblage script
	============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 23/05/2018
	:version: 0.1

	Script description
	------------------

	This program is used to check assembly quality

	Example
	-------

	>>> QualityAssemblage.py -d /homedir/user/work/dataAssemblt/ -o /homedir/user/work/result.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display QualityAssemblage.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the assembly files
						
		- \-o <path/to/output/file>, --outdirPath <path/to/output/file>
						path of the output file

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFasta , recupId,fasta2dict,sort_human


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''TThis program is used to check assembly quality ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display QualityAssemblage version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'path of directory that contains all the assembly files')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output file')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outFile= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)

if __name__ == "__main__":

	version="0.1" 

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in QualityAssemblage (Version " + version + ")       ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

################## Main ################################
	f = open(outFile,'w')
	f.write('n\tn:500\tL50\tmin\tN80\tN50\tN20\tE-size\tmax\tsum\tname\n')
	for files in os.listdir(directory):

		Pathfile = directory+files
		isFasta(Pathfile)
		strain = recupId(Pathfile.split('/')[-1])
		dico_fasta = fasta2dict(Pathfile)
		lengthGenome = 0
		nbScaffold = 0
		lengthN50 = 0
		lengthN80 = 0
		lengthN20 = 0
		Esize = 0
		L50 = 0
		n500 = 0
		first = True




		for elt in dico_fasta.values():
			nbScaffold += 1
			lengthGenome = lengthGenome +len(elt.seq)
			Esize = Esize + (len(elt.seq)*len(elt.seq))
			if len(elt.seq) > 500 :
				n500 += 1


		for elt  in sorted(dico_fasta.keys(), key=sort_human):
			if first :
				maxs =  len(dico_fasta[elt].seq)
				first = False 
			L50 += 1
			lengthN50 = lengthN50 + len(dico_fasta[elt].seq)
			if lengthN50 >= (lengthGenome/2) :
				break
		
		N50 = len(dico_fasta[elt].seq)

		for elt  in sorted(dico_fasta.keys(), key=sort_human):

			lengthN80 = lengthN80 + len(dico_fasta[elt].seq)
			if lengthN80 >= (lengthGenome*0.8) :
				break
		
		N80 = len(dico_fasta[elt].seq)

		for elt  in sorted(dico_fasta.keys(), key=sort_human):

			lengthN20 = lengthN20 + len(dico_fasta[elt].seq)
			if lengthN20 >= (lengthGenome*0.2) :
				break
		
		N20 = len(dico_fasta[elt].seq)
		Esize = int(Esize/lengthGenome)

		f.write('%s\t%s\t%s\t500\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(nbScaffold,n500,L50,N80,N50,N20,Esize,maxs,lengthGenome,strain))

	f.close()

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))


