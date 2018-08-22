#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display ABYSS_launch version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-r1', '--reads1',type = str, required=True, dest = 'r1', help = 'Path of  R1 file  which must be assembled')
	filesreq.add_argument('-r2', '--reads2',type = str, required=True, dest = 'r2', help = 'Path of R2 file  which must be assembled')
	filesreq.add_argument('-s', '--suffix',type = str, required=True, dest = 'suffix', help = 'Suffix of  R1 file  which must be assembled')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	r1 = os.path.abspath(args.r1)
	r2 = os.path.abspath(args.r2)
	suffix = args.suffix
	outDir= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	outDir = verifDir(outDir)
############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("        Welcome in ABYSS_launch (Version " + version + ")          ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

############# Main #########################

	isolate=r1.split('/')[-1 ].replace('%s'%suffix,'')
	print(form('\nLancement des jobs pour : '+isolate+'\n','green',['bold','underline']))
	for kmers in [20,30,40,50,60,70,80,90]:
		pathdirectory = outDir+'/abyss_assembly_'+isolate+'_'+str(kmers) + '/'
		os.system('mkdir -p '+pathdirectory)
		os.system("/usr/local/bioinfo/abyss/1.9.0/bin/abyss-pe name="+isolate+"_"+str(kmers)+" k="+str(kmers)+" in='"+r1+" "+r2+"' -o "+pathdirectory+"abyss_assembly_"+isolate+"_"+str(kmers)+".fasta")


############## end message ###########################²

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




