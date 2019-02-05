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
	-------

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
import vcf

# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human





def recupAC(id,record,MQmin,DPmin,QDmin):
    '''
    Function that allows to recover the nucleic acid according to the quality of mapping
    '''
    IUPAC = {"AT": "W", "CG": "S", "AC": "M", "GT": "K", "AG": "R", "CT": "Y","TA": "W", "GC": "S", "CA": "M", "TG": "K", "GA": "R", "TC": "Y"}

    if diploidy:
        for sample in record.samples:
            infoDiploidy = sample['GT']
    MQ = record.INFO['MQ']
    DP = record.INFO['DP']
    if DP < DPmin or MQ < MQmin:
        ac = 'N'*len(record.REF)
    elif record.ALT[0] == None :
        ac = record.REF
    elif record.ALT[0] != None and DP >= DPmin and MQ >= MQmin:
        QD =record.INFO['QD']
        if QD >= QDmin:
            if diploidy :
                if  infoDiploidy == '1/1' :
                    if len(record.REF) == 1 and len(record.ALT[0]) == 1:
                        ac = record.ALT[0]
                    elif len(record.REF) >= 2 :
                        ac = str(record.ALT[0]) + 'N'*(len(record.REF) -1)
                    elif len(record.ALT[0]) >= 2  :
                        ac = record.REF
                elif  infoDiploidy == '1/2' :
                    ac = IUPAC[str(record.ALT[0])+str(record.ALT[1])]
                elif infoDiploidy == '0/1' or '1/0' :
                    if len(record.REF) == 1 and len(record.ALT) == 1:
                        ac = IUPAC[str(record.REF)+str(record.ALT[0])]
                    elif len(record.REF) >= 2 :
                        ac = str(record.ALT[0]) + 'N'*(len(record.REF) -1)
                    elif len(record.ALT[0]) >= 2  :
                        ac = record.REF
            else :
                ac = record.ALT[0]
        else:
            ac = 'N'*len(record.REF)
            print("Warning, the gene {} have a SNP (QD = {}) which don't pass the filter (QD < {})".format(id,QD,QDmin))

    return str(ac)

if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to search a list of gene in a vcf file and give in output
                                     CDS,protein and gene fasta file.''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display vcfSearch version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-f', '--file', type=str, required=True, dest='vcf_file',
                          help='Path of vcf file which be used')

    filesreq.add_argument('-g', '--gff', type=str, required=False,default = 'None', dest='gff',
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
    vcf_file = os.path.abspath(args.vcf_file)
    listeGene = args.listeGene
    gff = args.gff
    prefix = args.prefix
    MQmin = args.MQmin
    QDmin = args.QDmin
    DPmin = args.DPmin
    diploidy = args.diploidy
    prcN = args.prcN
    ########### Gestion directory ##############
    verifFichier(vcf_file)
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
    if gff != 'None' :
        verifFichier(gff)
        gff = os.path.abspath(gff)
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
                            if position not in dicoRNA[Scaffold].keys() :
                                dicoRNA[Scaffold][position] = {}
                            if id not in dicoRNA[Scaffold].keys() :
                                dicoRNA[Scaffold][position][id] = [id,brin,'{0}:{1}-{2}'.format(Scaffold,start,end),start]

                    elif types == 'CDS' :
                        if id not in dicoCDS.keys():
                            dicoCDS[id] = []
                        dicoCDS[id].append([Scaffold,int(start),int(end)])
    else :
        print("\nYou choose to don't use gff file\n")
    nb = 0

    print('\nOpenning vcf file\n')
    dicoSeq = {}
    dicoInfo = {}
    listePosition = []
    # if gff != 'None':
    #     vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    #     for record in vcf_reader:
    #         if record.CHROM in dicoRNA.keys() and record.POS in dicoRNA[record.CHROM].keys():
    #             listeId = dicoRNA[record.CHROM][record.POS].keys()
    #             for id in listeId:
    #                 listePosition = [0] # ATTENTION PROBLEME DE listePOSITION réinitialisé comme avant
    #                 if record.POS in listePosition:
    #                     if len(record.REF) > len(ac):
    #                         ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
    #                         dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)
    #
    #                     else:
    #                         dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
    #                     for i in range(1, len(record.REF)):
    #                         listePosition.append(record.POS + i)
    #
    #                 # elif record.POS > int(listePosition[-1]) + 1 and id in dicoSeq.keys() :
    #                 #     dicoSeq[id] += 'N' * ((int(listePosition[-1]) + 1) - int(record.POS))
    #                 #
    #                 # elif id not in dicoSeq.keys() and int(dicoRNA[record.CHROM][record.POS][id][3]) < record.POS :
    #                 #     dicoSeq[id] = 'N' * (record.POS - int(record.POS - dicoRNA[record.CHROM][record.POS][id][3]))
    #
    #                 else:
    #                     listePosition = [record.POS]
    #                     brin = dicoRNA[record.CHROM][record.POS][id][1]
    #                     if record.INFO == {}:
    #                         ac = 'N'
    #                         if id not in dicoSeq.keys():
    #                             dicoSeq[id] = ac
    #                             dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
    #                         else:
    #                             dicoSeq[id] += ac
    #                     else:
    #                         brin = dicoRNA[record.CHROM][record.POS][id][1]
    #                         if id not in dicoSeq.keys():
    #                             ac = recupAC(id, record, MQmin, DPmin, QDmin)
    #                             dicoSeq[id] = ac
    #                             dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
    #                         else:
    #                             ac = recupAC(id, record, MQmin, DPmin, QDmin)
    #                             dicoSeq[id] += ac


    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.CHROM not in dicoSeq.keys() :
            if int(record.POS) != 1 :
                dicoSeq[record.CHROM] = 'N'*(record.POS- 1)

        elif old_chrom == record.CHROM and old_position +1 < int(record.POS) :
            dicoSeq[record.CHROM] += 'N'*(int(record.POS) - old_position)

        if record.INFO == {}:
            nb+= 1

        else :
            if record.CHROM not in dicoSeq.keys():
                dicoSeq[record.CHROM] = 'N' * nb
            else :
                dicoSeq[record.CHROM] += 'N'*nb
            nb = 0
            if record.POS in listePosition:
                if len(record.REF) > len(ac):
                    ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)
                else:
                    dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
                for i in range(1, len(record.REF)):
                    listePosition.append(record.POS + i)

            else:
                listePosition = [record.POS]
                id = record.CHROM
                if id not in dicoSeq.keys():
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] = ac
                else:
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] += ac

        old_position = int(record.POS)
        old_chrom = record.CHROM
    dicoSeq[record.CHROM] += 'N' * nb

    print('\nWrite fasta\n')
    if gff !='None' :
        with open('{0}_gene.fasta'.format(prefix),'w') as output_gene, \
                open('{0}_cds.fasta'.format(prefix), 'w') as output_cds,\
                open('{0}_protein.fasta'.format(prefix),'w') as output_prot, \
                open('{0}_protein_imcomplete.fasta'.format(prefix), 'w') as output_prot_imcomplete :
            for id in sorted(dicoSeq.keys() ,key = sort_human):
                sequence_chrom = dicoSeq[id]
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
                elif cds.count('N') <= (len(cds) * (prcN/100)):
                    prot = seqcds.translate()
                    descrip = '| position = {0}, length = {1}, brin = "{2}"'.format(dicoInfo[id][2], len(str(prot)), dicoInfo[id][1])
                    record = SeqRecord(prot, id=id + '_incomplete', name=id + '_incomplete', description=descrip)
                    SeqIO.write(record, output_prot_imcomplete, "fasta")

        with open('{0}_scaffold.fasta'.format(prefix),'w') as output_gene :

            for id in sorted(dicoSeq.keys() ,key = sort_human):
                seq = dicoSeq[id]
                record = SeqRecord(Seq(seq), id=id, name=id, description= f'| length : {len(seq)}')
                SeqIO.write(record, output_gene, "fasta")

