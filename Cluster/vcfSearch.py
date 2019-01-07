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
import argparse, os, sys,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from packages.progress.bar import ChargingBar
from collections import defaultdict

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human
from MODULES_SEB import AutoVivification




def recupAC(id,lineVCF,MQmin,DPmin,QDmin):
    '''
    Permet de récupéré l'acide nucléique à une position donnée en fonction de la qualité
    '''
    IUPAC = {"AT": "W", "CG": "S", "AC": "M", "GT": "K", "AG": "R", "CT": "Y","TA": "W", "GC": "S", "CA": "M", "TG": "K", "GA": "R", "TC": "Y"}
    if diploidy:
        vcf_scaffold, position, _, ref, alt, qual, _, info, *_,infoDiploidy = line.split()
        infoDiploidy = infoDiploidy.split(':')
        for elt in infoDiploidy:
            if '/' in elt :
                infoDiploidy = elt
                break
    else :
        vcf_scaffold, position, _, ref, alt, qual, _, info, *_ = line.split()

    if 'MQ=' in info and 'DP=' in info :
        MQ = float(info.split('MQ=')[1].split(';')[0])
        DP = float(info.split('DP=')[1].split(';')[0])
    else :
        MQ = 0
        DP = 0

    if DP < DPmin or MQ < MQmin:
        ac = 'N'
    elif alt == '.':
        ac = ref
    elif alt != '.' and DP >= DPmin and MQ >= MQmin:
        QD = float(info.split('QD=')[1].split(';')
        [0])
        if QD >= QDmin:
            if diploidy :
                if  infoDiploidy == '1/1' :
                    ac = alt
                elif  infoDiploidy == '1/2' :
                    ac = IUPAC[alt.replace(',','')]
                elif infoDiploidy == '0/1' or '1/0' :
                    if len(ref) == 1 and len(alt) == 1:
                        ac = IUPAC[ref+alt]
                    elif len(ref) >= 2 :
                        ac = alt + 'N'*(len(ref) -1)
                    elif len(alt) >= 2  :
                        ac = ref
            else :
                ac = alt
        else:
            ac = 'N'*len(ref)
            print("Warning, the gene {} have a SNP (QD = {}) which don't pass the filter (QD < {})".format(id,QD,QDmin))

    return ac

if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to search a list of gene in a vcf file and give in output
                                     CDS,protein and gene fasta file.''')
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
    files.add_argument('-MQ', '--MQ', type=int, required=False,default = 20, dest='MQmin', help='The minimum mapping quality supporting to accept the mapping (default = 20)')
    files.add_argument('-DP', '--DP', type=int, required=False, default= 10, dest='DPmin',
                          help=' The minimum sequencing depth supporting to accept the mapping ( default = 10)')
    files.add_argument('-QD', '--QD', type=int, required=False, default= 5, dest='QDmin',
                          help='The minimum QD value to accept the variant ( default = 5)')
    files.add_argument('-prop', '--proportion', type=int, required=False, default= 80, dest='prcN',
                       help='The maximum percentage of N in a sequence (default = 20)')
    files.add_argument('-diploidy', '--diploidy',action='store_true', dest='diploidy',
                          help='Use this option if you are using a diploid organism')
    files.add_argument('-p', '--prefix', type=str, required=False, default='vcfSearch', dest='prefix',
                          help='prefix for the output file')

    ######### Recuperation arguments ###########
    args = parser.parse_args()
    vcf = os.path.abspath(args.vcf)
    listeGene = args.listeGene
    gff = os.path.abspath(args.gff)
    prefix = args.prefix
    MQmin = args.MQmin
    QDmin = args.QDmin
    DPmin = args.DPmin
    diploidy = args.diploidy
    prcN = args.prcN
    ########### Gestion directory ##############
    verifFichier(vcf)
    verifFichier(gff)
    if listeGene != 'None' :
        verifFichier(listeGene)
        listeGene = os.path.abspath(listeGene)
    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
        "\t" + form("|", 'yellow', 'bold') + form("            Welcome in vcfSearch (Version " + version + ")         ",
                                                  type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

    ########### Main #####################################
    liste = []
    if listeGene != 'None' :
        lines = openfile(listeGene)
        liste = [elt.split()[0] for elt in lines]

    dicoRNA = {}
    dicoCDS = {}
    print('\nOpenning gff file\n')

    with open(gff,'r') as gff_file :
        header = gff_file.readline()
        for line in gff_file :
            try :
                Scaffold, _ ,types,start,end,_ ,brin,_,id = line.rstrip().split('\t')
            except :
                print(line)
                exit()
            id = id.split(';')[0].split('=')[1]
            if ':' in id :
                id = id.split(':')[0]
            if id in liste or listeGene == 'None' :
                if types == 'mRNA' :
                    for position in range(int(start),int(end)+1) :
                        if Scaffold not in dicoRNA.keys():
                            dicoRNA[Scaffold] = {}
                        if id not in dicoRNA[Scaffold].keys() :
                            dicoRNA[Scaffold][position] = [id,brin,'{0}:{1}-{2}'.format(Scaffold,start,end),start]
                elif types == 'CDS' :
                    if id not in dicoCDS.keys():
                        dicoCDS[id] = []
                    dicoCDS[id].append([Scaffold,int(start),int(end)])

    print('\nOpenning vcf file\n')
    dicoSeq = {}
    dicoInfo = {}
    listePosition = []
    with open(vcf, "r") as f:
        for line in f :
            if line[0] != '#' :
                vcf_scaffold, position, _, ref, alt, qual, _, info, *_, = line.split()
                position = int(position)
                if vcf_scaffold in dicoRNA.keys() and position in dicoRNA[vcf_scaffold].keys():
                    if position in listePosition :
                        dicoSeq[id] = dicoSeq[id].replace(dicoSeq[id][-len(ac):],'N'*len(ac))
                        for i in range(1,len(ref)) :
                            listePosition.append(position+i)

                    else :
                        listePosition = [position]
                        id = dicoRNA[vcf_scaffold][position][0]
                        brin = dicoRNA[vcf_scaffold][position][1]
                        if id not in dicoSeq.keys():
                            ac = recupAC(id, line, MQmin, DPmin, QDmin)
                            dicoSeq[id] = ac
                            dicoInfo[id] = dicoRNA[vcf_scaffold][position]
                        else :
                            ac = recupAC(id, line, MQmin, DPmin, QDmin)
                            dicoSeq[id] += ac



    print('\nWrite fasta\n')
    with open('{0}_gene.fasta'.format(prefix),'w') as output_gene, \
            open('{0}_cds.fasta'.format(prefix), 'w') as output_cds,\
            open('{0}_protein.fasta'.format(prefix),'w') as output_prot, \
            open('{0}_protein_imcomplete.fasta'.format(prefix), 'w') as output_prot_imcomplete :

        for id in sorted(dicoSeq.keys() ,key = sort_human):
            seq = dicoSeq[id]
            brin =  dicoInfo[id][1]
            start  = dicoInfo[id][3]
            descrip = '| position = {0}, length = {1}, brin = "{2}"'.format(dicoInfo[id][2], len(str(seq)), dicoInfo[id][1])
            record = SeqRecord(Seq(seq), id=id, name=id, description=descrip)
            SeqIO.write(record, output_gene, "fasta")
            cds = ''
            for elt in dicoCDS[id]:
                cds = cds + seq[elt[1] - (int(start)):(elt[2] +1 - int(start))]
            if brin == '-':
                seqcds = Seq(cds).reverse_complement()
            elif brin == '+':
                seqcds = Seq(cds)
            descrip = '| position = {0}, length = {1}, brin = "{2}"'.format(dicoInfo[id][2], len(str(cds)), dicoInfo[id][1])
            record = SeqRecord(seqcds, id=id, name=id, description=descrip)
            SeqIO.write(record, output_cds, "fasta")
            if 'N' not in cds:
                prot = seqcds.translate()
                descrip = '| position = {0}, length = {1}, brin = "{2}"'.format(dicoInfo[id][2], len(str(prot)), dicoInfo[id][1])
                record = SeqRecord(prot, id=id, name=id, description=descrip)
                SeqIO.write(record, output_prot, "fasta")
            elif cds.count('N') < (len(cds) * (prcN/100)):
                prot = seqcds.translate()
                descrip = '| position = {0}, length = {1}, brin = "{2}"'.format(dicoInfo[id][2], len(str(prot)), dicoInfo[id][1])
                record = SeqRecord(prot, id=id + '_incomplete', name=id + '_incomplete', description=descrip)
                SeqIO.write(record, output_prot_imcomplete, "fasta")


