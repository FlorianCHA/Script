#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package blast2fasta.py
# @author Florian CHARRIAT

"""
    The blast2fasta script
    ===========================
    :author: Florian CHARRIAT
    :contact: florian.charriat@inra.fr
    :date: 08/07/2018
    :version: 0.1

    Script description
    ------------------
    This Programme is used to take the blast result and give a fasta format with the environnement of the sequence

    Example
    -------
    >>> blast2fasta.py -b blast_resut.txt -f assembly.fasta -l 1500:200 -o result.fasta


    Help Programm
    -------------
    optional arguments:
        - \-h, --help
                        show this help message and exit
        - \-v, --version
                        display blast2fasta.py version number and exit
    Input mandatory infos for running:
        - \-b <path/to/blast/file>, --blast <path/to/directory>
                        path of result blast result with the option outfmt 6
        - \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
                        path of fasta file which be used for the blast (fasta assembly)
        - \-l <str>, --len <str>
                        the size upstream and downstream of the sequence (at format ustream:downstream) by default both values ​​are 0.
        - \-o  <path/to/output/fasta/file>, --len <path/to/output/fasta/file>
                        path of the output fasta file

"""


##################################################
## Modules
##################################################
#Import
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human, re

## Python modules
import argparse

## BioPython modules
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
    version =  '0.1'

    # Parameters recovery
    parser = argparse.ArgumentParser(prog=__file__, description='''This Programme i used to search ET in gene promoter''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
                        'display '+__file__+' version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-b', '--blast',type = str, required=True, dest = 'blast', help = 'path of result blast result with the option outfmt 6')
    filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = ' path of fasta file which be used for the blast (fasta assembly)')
    filesreq.add_argument('-l', '--length', type=str, required=False, default= '0:0', dest='length', help='the size upstream and downstream of the sequence (at format ustream:downstream) by default both values ​​are 0')
    filesreq.add_argument('-o', '--output', type=str, required=True ,dest='output', help='path of the output fasta file')


######### Recuperation arguments ###########
    args = parser.parse_args()
    fasta = os.path.abspath(args.fasta)
    blast = os.path.abspath(args.blast)
    length = args.length
    output = os.path.abspath(args.output)
    verifFichier(fasta)
    verifFichier(blast)
    up, down = length.split(':')
    if up == '' :
        up = 0
    else :
        up = int(up)
    if down == '' :
        down = 0
    else :
        down = int(down)
########### Gestion directory ##############

    ############### start message #########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
            "\t" + form("|", 'yellow', 'bold') + form(
        "            Welcome in blast2fasta (Version " + version + ")         ",
        type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')
    
    ################ Main ################################

    ## A dico is create with all sequence of the fasta file ##
    print(form(f'{fasta.split("/")[-1]} file processing\n','green','bold'))
    fasta_dico = fasta2dict(fasta)

    ## The blast file is read and the output is create

    print(form(f'{blast.split("/")[0]} file processing & Write the output file\n'))
    with open(blast,'r') as blast_file, open(output,'w') as output_file :
        for line in blast_file :
            start = int(line.split()[8])
            stop = int(line.split()[9])
            scaffold = line.split()[1]
            if start > stop :
                seq = fasta_dico[scaffold].seq
                seq = Seq(str(seq)[(stop -1  - down):(start + up)])
                seq = seq.reverse_complement()
            else :
                seq = Seq(str(fasta_dico[scaffold].seq)[(start -1 - up):(stop+down)])

            record = SeqRecord(seq, id=str(fasta.split("/")[-1]), name=str(fasta.split("/")[-1]), description='| length : ' + str(len(seq)))
            SeqIO.write(record, output_file, "fasta")


    

    ############## end message ###########################

    print(form("\n\t---------------------------------------------------------",'yellow','bold'))
    print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
    print(form("\t---------------------------------------------------------",'yellow','bold'))