#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package triFasta.py
# @author Florian CHARRIAT

"""
    The triFasta script
    ===========================
    :author: Florian CHARRIAT
    :contact: florian.charriat@inra.fr
    :date: 08/07/2018
    :version: 0.1

    Script description
    ------------------
    This program is used to delete a sequence of less than or greater than a given length and a sequence with more than x N, where x is given by the user

    Example
    -------
    >>> triFasta.py -b blast_resut.txt -f assembly.fasta -l 1500:200 -o result.fasta


    Help Programm
    -------------
    optional arguments:
        - \-h, --help
                        show this help message and exit
        - \-v, --version
                        display triFasta.py version number and exit
    Input mandatory infos for running:
        - \-f <path/to/fasta/file>, --fasta <path/to/fasta/file>
                        path of fasta file which be used.
        - \-l <int>, --len <int>
                        Lensize cutoff ( use only if you use the k option)
        - \-k <g/greater/l/lower>, --keep <g/greater/l/lower>
                        Choice keep sequences size greater than -l (g/greater) or keep lower (l/lower)
        - \-n <int>, --numberN int>
                        Number max of N nuclétotide in sequence
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
    filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = ' path of fasta file which be used')
    filesreq.add_argument('-l', '--len', type=int, required=False, default= 0, dest='length', help='Lensize cutoff ( use only if you use the k option)')
    filesreq.add_argument('-k', '--keep', type=str, required=False, default='None', dest='keep',
                          help='Choice keep sequences size greater than -l (g/greater) or keep lower (l/lower)')
    filesreq.add_argument('-n', '--numberN', type=int, required=False, default= 0, dest='numberN',
                          help='Number max of N nuclétotide in sequence')
    filesreq.add_argument('-o', '--output', type=str, required=True ,dest='output', help='path of the output fasta file')


######### Recuperation arguments ###########
    args = parser.parse_args()
    fasta = os.path.abspath(args.fasta)
    length = args.length
    k = args.keep
    numberN = args.numberN
    output = os.path.abspath(args.output)
    verifFichier(fasta)
########### Gestion directory ##############

    ############### start message #########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
            "\t" + form("|", 'yellow', 'bold') + form(
        "            Welcome in triFasta (Version " + version + ")         ",
        type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')
    
    ################ Main ################################

    ## A dico is create with all sequence of the fasta file ##
    print(form(f'{fasta.split("/")[-1]} file processing\n','green','bold'))
    fasta_dico = fasta2dict(fasta)
    with open(output,'w') as output_file :
        for elt in sorted(fasta_dico.keys(), key= sort_human) :

            if str(fasta_dico[elt].seq).count('N') > numberN:
                print(f'{elt} has removed, they havec {str(fasta_dico[elt].seq).count("N")} N in sequence ! \n')

            elif k == 'g' or k == 'greater' :
                if len(fasta_dico[elt].seq) > length :
                    record = SeqRecord(fasta_dico[elt].seq, id=str(elt), name=str(elt),
                                       description='| length : ' + str(len(fasta_dico[elt].seq)))
                    SeqIO.write(record, output_file, "fasta")
                else :
                    print(f'{elt} has removed, the length is lower than {length}! \n')

            elif k == 'l' or k == 'lower':
                if len(fasta_dico[elt].seq) < length:
                    record = SeqRecord(fasta_dico[elt].seq, id=str(elt), name=str(elt),
                                       description='| length : ' + str(len(fasta_dico[elt].seq)))
                    SeqIO.write(record, output_file, "fasta")
                else:
                    print(f'{elt} has removed ! \n')

            else :
                record = SeqRecord(fasta_dico[elt].seq, id=str(elt), name=str(elt),
                                   description='| length : ' + str(len(fasta_dico[elt].seq)))
                SeqIO.write(record, output_file, "fasta")






    

    ############## end message ###########################

    print(form("\n\t---------------------------------------------------------",'yellow','bold'))
    print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
    print(form("\t---------------------------------------------------------",'yellow','bold'))