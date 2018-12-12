#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package vcfSearch.py
# @author Florian Charriat

"""
	The vcfSearch script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 13/11/2018
	:version: 0.1

	Script description
	------------------

	This program is used to search a list of gene in a vcf file and give in output CDS,protein and gene fasta file.

	Example
	-------qstat

	>>> vcfSearch.py -l /homedir/user/work/listegene.txt -g /homedir/user/work/Isolat.gff -f /homedir/user/work/file.vcf -p prefix

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display vcfSearch.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/vcf/file>, --file <path/to/vcf/file>
						path of vcf file which be used

		- \-l <path/to/gene/liste/file>, --listegene <path/to/gene/liste/file>
						path of the file which contains a gene list of interest. If you don't give a list,
						the script search all gene of the gff file

        - \-g <path/to/gff/file>, --gff <path/to/gff/file>
						path of the gff file path of the reference genome used to create the vcf

		- \-p <prefix>, --prefix <prefix>
						prefix for the output file

"""

########## Module ###############
## Python modules
import argparse, os, sys, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from packages.progress.bar import ChargingBar
from collections import defaultdict

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human
from MODULES_SEB import AutoVivification



def recupstat(file,line,map,noMap):
    '''
    Permet de récupéré l'acide nucléique à une position donnée en fonction de la qualité
    '''
    vcf_scaffold, position, _, ref, alt, qual, _, info, *_ = line.split()

    if 'MQ=' in info and 'DP=' in info:
        MQ = float(info.split('MQ=')[1].split(';')[0])
        DP = float(info.split('DP=')[1].split(';')[0])
        if alt == '.' :
            QD = 0
        elif "QD=" in line :
            QD = float(info.split('QD=')[1].split(';')[0])
        file.write(f'{position}\t{MQ}\t{DP}\t{QD}\n')
        map += 1
    else:
        noMap += 1

    return map,noMap


if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to get mapping statistcs from vcf file of  list of gene''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display vcfSearch version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-f', '--file', type=str, required=True, dest='vcf',
                          help='Path of vcf file which be used')

    filesreq.add_argument('-g', '--gff', type=str, required=True, dest='gff',
                          help='Path of the gff file path of the reference genome used to create the vcf')

    files = parser.add_argument_group('Input infos for running with default values')
    files.add_argument('-l', '--listegene', type=str, required=False, default='None', dest='listeGene',
                       help='Path of the file which contains a gene list of interest,'
                            'the script search all gene of the gff file')
    files.add_argument('--cds',action='store_true', dest='cds',
                       help='Retrieve only the stat of CDS mapping')
    files.add_argument('--gene',action='store_true', dest='gene',
                       help='Retrieve only the stat of gene mapping')
    files.add_argument('-p', '--prefix', type=str, required=False, default='vcfSearch', dest='prefix',
                       help='prefix for the output file')

    ######### Recuperation arguments ###########
    args = parser.parse_args()
    vcf = os.path.abspath(args.vcf)
    listeGene = args.listeGene
    gff = os.path.abspath(args.gff)
    prefix = args.prefix
    cds = args.cds
    gene = args.gene

    if gene == False and cds == False or gene == True and cds == True :
        print(form("\nWarning you must choose between gene option and cds option.\n",'orange','bold'))
        exit()

    ########### Gestion directory ##############
    verifFichier(vcf)
    verifFichier(gff)
    if listeGene != 'None':
        verifFichier(listeGene)
        listeGene = os.path.abspath(listeGene)
    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
            "\t" + form("|", 'yellow', 'bold') + form(
        "            Welcome in vcfSearch (Version " + version + ")         ",
        type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

    ########### Main #####################################
    liste = []
    if listeGene != 'None':
        lines = openfile(listeGene)
        liste = [elt.split()[0] for elt in lines]

    dicoRNA = {}
    dicoCDS = {}
    print('\nOpenning gff file\n')

    with open(gff, 'r') as gff_file:
        header = gff_file.readline()
        for line in gff_file:
            try:
                Scaffold, _, types, start, end, _, brin, _, id = line.rstrip().split('\t')
            except:
                print(line)
                exit()
            id = id.split(';')[0].split('=')[1]
            if ':' in id:
                id = id.split(':')[0]
            if id in liste or listeGene == 'None':
                if types == 'mRNA':
                    if gene == True :
                        for position in range(int(start), int(end) + 1):
                            if Scaffold not in dicoRNA.keys():
                                dicoRNA[Scaffold] = {}
                            if id not in dicoRNA[Scaffold].keys():
                                dicoRNA[Scaffold][position] = [id, brin, '{0}:{1}-{2}'.format(Scaffold, start, end), start]
                elif types == 'CDS':
                    if cds == True :
                        for position in range(int(start), int(end) + 1):
                            if Scaffold not in dicoRNA.keys():
                                dicoRNA[Scaffold] = {}
                            if id not in dicoRNA[Scaffold].keys():
                                dicoRNA[Scaffold][position] = [id, brin, '{0}:{1}-{2}'.format(Scaffold, start, end), start]

    print('\nOpenning vcf file\n')
    dicoSeq = {}
    dicoInfo = {}
    map = noMap = 0
    with open(vcf, "r") as f,open(f'{prefix}_value.csv','w') as value_file:
        value_file.write('Position\tMQ\tDP\tQD\n')
        for line in f:
            if line[0] != '#':
                vcf_scaffold, position, _, ref, alt, qual, _, info, *_, = line.split()
                position = int(position)
                if vcf_scaffold in dicoRNA.keys() and position in dicoRNA[vcf_scaffold].keys():
                    id = dicoRNA[vcf_scaffold][position][0]
                    brin = dicoRNA[vcf_scaffold][position][1]
                    map,noMap = recupstat(value_file,line,map,noMap)


    print(form(f'\n\tThere is {map} position mapped to {map+noMap} position ({round((map/(map+noMap))*100,2)}%)\n','green','bold'))
