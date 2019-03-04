#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Merge_fasta.py
# @author Florian Charriat

"""
	The Merge_fasta script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 13/11/2018
	:version: 0.1

	Script description
	------------------

	This program is used to merge fasta file and sort Sequnce ID

	Example
	-------

	>>> Merge_fasta.py -f file.fasta -o output.fasta

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display Merge_fasta.py version number and exit

	Input mandatory infos for running:
		- --fasta1 <path/to/fasta/file>
						path of the first fasta file which be used

		- \--fasta2 <path/to/fasta/file>
						path of the second fasta file which be used

		- \-o  <path/to/output/file>, --output  <path/to/output/file>
						path of the output file

"""

########## Module ###############
## Python modules
import argparse, os, sys,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human

if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to merge fasta file and sort Sequnce ID.''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display gvcfSearch version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('--fasta1', type=str, required=True, dest='fasta1',
                          help='Path of the first fasta file which be used')
    filesreq.add_argument('--fasta2', type=str, required=True, dest='fasta2',
                          help='Path of the second fasta file which be used')
    filesreq.add_argument('-o', '--output', type=str, required=True, dest='output',
                          help='Path of the output file')

    ######### Recuperation arguments ###########
    args = parser.parse_args()
    fasta1 = os.path.abspath(args.fasta1)
    fasta2 = os.path.abspath(args.fasta2)
    output = os.path.abspath(args.output)
    ########### Gestion directory ##############
    verifFichier(fasta1)
    verifFichier(fasta2)

    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
            "\t" + form("|", 'yellow', 'bold') + form(
        "            Welcome in alnToSnp (Version " + version + ")         ",
        type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

########### Main #####################################
    seq1 = fasta2dict(fasta1)
    seq2 = fasta2dict(fasta2)
    liste = list(seq1.keys())
    liste = liste + (list(seq2.keys()))
    liste = list(set(liste))
    with open(output, 'w') as output_file :
        for elt in sorted(liste, key = sort_human) :
            if elt in seq1.keys():
                sequence = seq1[elt].seq
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description= seq1[elt].description)
                SeqIO.write(record, output_file, "fasta")
            if elt in seq2.keys():
                sequence = seq2[elt].seq
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description= seq2[elt].description)
                SeqIO.write(record, output_file, "fasta")