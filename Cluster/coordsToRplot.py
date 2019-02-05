#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package coordsToRplot.py
# @author Florian Charriat

"""
	The coordsToRplot script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 13/11/2018
	:version: 0.1

	Script description
	------------------

	This program is used to convert show-coords output (MUMmer tools) to file for genoPlotR package

	Example
	-------

	>>> coordsToRplot.py -l liste_scaffold.txt -f show_coords.output -p output.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display coordsToRplot.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/show_coords/file>, --file <path/to/show_coords/file>
						path of show_coords file which be used

		- \-l <path/to/scaffolds/liste/file>, --liste <path/to/scaffolds/liste/file>
						path of the file which contains a scaffold list with length ( Scaffold : 4000000)

		- \-p <prefix>, --prefix <prefix>
						prefix for the output file

"""

########## Module ###############
## Python modules
import argparse, os, sys,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import vcf

from module_Flo import fasta2dict,verifFichier,form,sort_human

if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to search a list of gene in a vcf file and give in output
                                     CDS,protein and gene fasta file.''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display coordsToRplot version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-f', '--file', type=str, required=True, dest='file',
                          help='path of show_coords file which be used')
    filesreq.add_argument('--fasta1', type=str, required=False, default='None', dest='fasta1',
                          help='path of the ref assembly used with MUMmer')
    filesreq.add_argument('--fasta2', type=str, required=False, default='None', dest='fasta2',
                          help='path of the query assembly used with MUMmer')
    filesreq.add_argument('-p', '--prefix', type=str, required=False, default='prefix', dest='prefix',
                          help='prefix for the output file')

    ######### Recuperation arguments ###########
    args = parser.parse_args()
    file = os.path.abspath(args.file)
    fasta1 = os.path.abspath(args.fasta1)
    fasta2 = os.path.abspath(args.fasta2)
    prefix = os.path.abspath(args.prefix)
    ########### Gestion directory ##############
    verifFichier(file)
    verifFichier(fasta1)
    verifFichier(fasta2)
    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
        "\t" + form("|", 'yellow', 'bold') + form("            Welcome in coordsToRplot (Version " + version + ")         ",
                                                  type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

    ########### Main #####################################
    print('\nChargement des fasta (1/2)')
    f1 = fasta2dict(fasta1)
    print('\nChargement des fasta (2/2)')
    f2 = fasta2dict(fasta2)
    dico_len = {}
    for elt in f1.keys() :
        dico_len[elt] = len(f1[elt].seq)
    for elt in f2.keys() :
        dico_len[elt] = len(f2[elt].seq)
    print('\nFinish')
    dico = {}
    with open(file,'r') as input_file :
        for line in input_file :
            s1,e1,_,s2,e2,*_,c1,c2 = line.split()
            if c1 not in dico.keys():
                dico[c1] = [[int(s1),int(e1),int(s2),int(e2),c2]]
            else :
                dico[c1].append([int(s1),int(e1),int(s2),int(e2),c2])
    ajout1 = ajout2 = 0
    old_c1 = sorted(dico.keys(), key = sort_human)[0]
    old_c2 = dico[sorted(dico.keys(), key = sort_human)[0]][0][-1]
    dico1 = {}
    with open(f'{prefix}_dnaseg_ref.txt','w') as dnaR_file , open(f'{prefix}_limits_ref.txt','w') as limitR_file:
        for elt in sorted(dico.keys(), key = sort_human):
            liste = sorted(dico[elt])
            for coords in liste :
                s1, e1, s2, e2, c2 = coords
                if old_c1 != elt:
                    dnaR_file.write(f'{old_c1}\t{ajout1 + 1}\t{ajout1 + dico_len[old_c1]}\t1\n')
                    ajout1 = ajout1 + dico_len[old_c1]
                    limitR_file.write(f'{ajout1}\n')
                if c2 not in dico1.keys() :
                    dico1[c2] = [[s1+ajout1,e1 +ajout1,s2,e2,c1]]
                else :
                    dico1[c2].append([s1 + ajout1, e1 + ajout1, s2, e2,c1])


                old_c1 = elt
        dnaR_file.write(f'{old_c1}\t{ajout1 + 1}\t{ajout1 + dico_len[old_c1]}\t1\n')

    with open(f'{prefix}_comparaison.txt','w') as comparaison_file, open(f'{prefix}_limits_query.txt','w') as limitQ_file,open(f'{prefix}_dnaseg_query.txt','w') as dnaQ_file :
        for elt in sorted(dico1.keys(), key=sort_human):
            liste = sorted(dico1[elt])
            for coords in liste:
                s1, e1, s2, e2, c1 = coords
                if old_c2 != elt:
                    dnaQ_file.write(f'{old_c2}\t{ajout2 + 1}\t{ajout2 + dico_len[old_c2]}\t1\n')
                    ajout2 = ajout2 + dico_len[old_c2]
                    limitQ_file.write(f'{ajout2}\n')
                old_c2 = elt
                comparaison_file.write(f'{s1}\t{e1}\t{s2 + ajout2}\t{e2 + ajout2}\n')


        dnaQ_file.write(f'{old_c2}\t{ajout2 + 1}\t{ajout2 + dico_len[old_c2]}\t1\n')
