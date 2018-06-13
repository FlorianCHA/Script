#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package renameAndAddParent.py
# @author Sebastien Ravel



##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import dict2txt, loadInList, dictDict2txt, existant_file,fasta2dict, sort_human

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='26/06/2017'


##################################################

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme rename genes name''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input  infos for running')
	filesreq.add_argument('-g', '--gff', metavar="<path/to/gff>", type=existant_file, required=True, dest = 'gffFileIn', help = 'path to gff file')
	filesreq.add_argument('-s', '--strain', metavar="<str>", type=str, required=True, dest = 'strainName', help = 'Name of strain')
	filesreq.add_argument('-o', '--out', metavar="<path/to/output/file>", required=True, dest = 'outDir', help = 'Path of output file')
	filesreq.add_argument('-n', '--numero', metavar="<int>", type=int, required=False,default=0 ,dest = 'num', help = 'path to gff file')
	filesreq.add_argument('-t', '--tools', metavar="<str>", type=str, required=False,default='Braker', dest = 'tools', help = 'path to gff file')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome in renameGFF (Version " + version + ")               #")
	print("#################################################################\n\n")


	gffFileIn = args.gffFileIn
	strainName = args.strainName
	ouputDir = args.outDir
	num =args.num
	tools = args.tools
	# resume value to user
	print(" - Intput Info:")
	print("\t - GFF file is : %s" % gffFileIn)
	print("\t - Strain Name is : %s" % strainName)

	#print(" - Output Info:")
	#print("\t - Output fasta with align is:  %s\n\n" % outFile)

	# parse le GFF
	#Chromosome_8.1	AUGUSTUS	gene	1	4965	0.1	-	.	ID=g1;
	#Chromosome_8.1	AUGUSTUS	mRNA	1	4965	0.1	-	.	ID=g1.t1;Parent=g1
	#Chromosome_8.1	AUGUSTUS	CDS	1	911	0.43	-	2	ID=g1.t1.CDS1;Parent=g1.t1
	#Chromosome_8.1	AUGUSTUS	exon	1	911	.	-	.	ID=g1.t1.exon1;Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	intron	912	972	0.44	-	.	Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	CDS	973	2883	0.19	-	2	ID=g1.t1.CDS2;Parent=g1.t1
	#Chromosome_8.1	AUGUSTUS	exon	973	2883	.	-	.	ID=g1.t1.exon2;Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	intron	2884	2994	0.43	-	.	Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	CDS	2995	3461	0.28	-	1	ID=g1.t1.CDS3;Parent=g1.t1
	#Chromosome_8.1	AUGUSTUS	exon	2995	3461	.	-	.	ID=g1.t1.exon3;Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	intron	3462	3631	0.59	-	.	Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	CDS	3632	4965	0.6	-	0	ID=g1.t1.CDS4;Parent=g1.t1
	#Chromosome_8.1	AUGUSTUS	exon	3632	4965	.	-	.	ID=g1.t1.exon4;Parent=g1.t1;
	#Chromosome_8.1	AUGUSTUS	start_codon	4963	4965	.	-	0	Parent=g1.t1;
	startPoly=0
	with open(gffFileIn, "r") as gffFile, open(ouputDir,"w") as outFileGff:
		if num ==0 :
			outFileGff.write("##gff-version 3\n")
		for line in gffFile:
			if line[0] != "#":

				tabLine = line.rstrip().split("\t")
				typeSeq = tabLine[2]
				if tabLine[1] == "AUGUSTUS":

					# GENE
					if typeSeq == "gene":
						firstCDS = 0
						# new gene so write polyInfo
						if startPoly != 0:
							polyLine = "{0}\tpolypeptide\t{1}\t{2}\t{3}\tID={4};Derives_from={5};Name={4}\n".format("\t".join(tabLine[:2]),startPoly, stopPoly, "\t".join(tabLine[5:8]), polyName, mRNAName)
							outFileGff.write(polyLine.replace("AUGUSTUS",tools+'_BGPI'))

						geneNumero = tabLine[8].split(";")[0].replace("ID=","").replace("g","")
						geneNumeroReformat = geneNumero
						if num != 0 :
							num = num +1
							geneNumeroReformat = str(num)
						geneNumeroReformat = str(geneNumeroReformat.zfill(5)) +'0'
						geneName = tabLine[8].split(";")[0].replace("g"+geneNumero,"Mo_"+strainName+"_"+geneNumeroReformat).replace("ID=","")

						#print(geneNumero)
						#print(geneNumeroReformat)
						#print(geneName)

						newline = "{0}\tID={1};Name={1}\n".format("\t".join(tabLine[:8]), geneName)
						outFileGff.write(newline.replace("AUGUSTUS",tools+'_BGPI'))

					# transcript => mRNA
					elif typeSeq == "mRNA" or typeSeq == 'transcript':
						numT = str(tabLine[8].split(";")[0].split('.t')[-1])
						NewnumT = str(int(numT)-1)
						mRNAName = tabLine[8].split(";")[0].replace("g"+geneNumero,"Mo_"+strainName+"_"+geneNumeroReformat).replace("ID=","").replace('.t%s'%numT,'T%s'%NewnumT)
						newline = "{0}\tID={1};Parent={2}\n".format("\t".join(tabLine[:8]), mRNAName, geneName)
						outFileGff.write(newline.replace("transcript","mRNA").replace("AUGUSTUS",tools+'_BGPI'))

					# CDS
					elif typeSeq == "CDS" :
						stopPoly = tabLine[4]
						CDSLine = "%s\tParent=%s\n" % ("\t".join(tabLine[:8]), mRNAName)
						outFileGff.write(CDSLine.replace("AUGUSTUS",tools+'_BGPI'))

					# exon ID=g1.t1.exon2;Parent=g1.t1;
					elif typeSeq == "exon":
						exonLine = "%s\tParent=%s\n" % ("\t".join(tabLine[:8]),mRNAName)
						outFileGff.write(exonLine.replace("AUGUSTUS",tools+'_BGPI'))

					# intron start_codon stop_codon
					elif typeSeq == "intron":
						intronLine = "%s\tParent=%s\n" % ("\t".join(tabLine[:8]), mRNAName)
						outFileGff.write(intronLine.replace("AUGUSTUS",tools+'_BGPI'))

					elif typeSeq == "start_codon" or typeSeq == "stop_codon":
						Line = "%s\tParent=%s\n" % ("\t".join(tabLine[:8]), mRNAName)
						outFileGff.write(Line.replace("AUGUSTUS",tools+'_BGPI'))

				else:
					outFileGff.write(line.replace("AUGUSTUS",tools+'_BGPI'))
			else:
				outFile.write(line.replace("AUGUSTUS",tools+'_BGPI'))

