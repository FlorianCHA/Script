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


########## Fonction ###############
def sort_human(s, _nsre=re.compile('([0-9]+)')):
    """ Sort the list in the way that humans expect, use list.sort(key=sort_human) or sorted(list, key=sort_human)).
    """
    try:
        return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]
    except TypeError:
        if not isinstance(s, int):
            print("WARNNING MODULES_SEB::sort_human : List %s value not understand so don't sort \n" % s)


    return s

def verifFichier(fichier):
    '''Permet de vérifier si un fichier existe.

    :Parameters:
         fichier
        Path du fichier

    '''
    if os.path.exists(fichier):
        return
    else:
        raise ValueError(form("ERROR the file '%s' doesn't exist, please check if your files exists" % fichier, 'red', 'bold'))

def openfile(file):
    """
	Permet de mettre dans une varaible une liste contenant chaque ligne du fichier
	:Parameters:
		file
		  Le chemin du fichier à ouvrir
	"""
    f = open(file,'r')
    lines = f.readlines()
    f.close()
    return lines

def fasta2dict(filename):
    """
	Function that take a file name (fasta), and return a dictionnary of sequence
	"""
    with open(filename, "rU") as fastaFile:
        return SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))

def form(text, col='white', type='none'):
    '''
    Permet de mettre en forme les textes afficher sur le terminale.

    :Parameters:
         text
        Le texte à transformer
         col
        La couleur souhaité entre les couleurs red, green, yellow, orange, blue et purple
         text
         str ou liste de str du format à appliquer (bold, underline, blind et highligth)
    '''
    W = '\033[0'  # white (normal)
    R = '\033[31'  # red
    G = '\033[32'  # green
    Y = '\033[33'  # yellow
    O = '\033[33'  # orange
    B = '\033[34'  # blue
    P = '\033[35'  # purple
    end = '\033[0m'  # white (normal)
    Bold = ';1'
    underline = ';4'
    blind = ';5'
    highlight = ';7'
    text = 'm' + text
    if 'bold' in type:
        text = Bold + text
    if 'underline' in type:
        text = underline + text
    if 'highlight' in type:
        text = blind + text
    if 'highlight' in type:
        text = highlight + text
    if col == 'red':
        return R + text + end
    elif col == 'white':
        return W + text + end
    elif col == 'green':
        return G + text + end
    elif col == 'yellow':
        return Y + text + end
    elif col == 'orange':
        return O + text + end
    elif col == 'blue':
        return B + text + end
    elif col == 'purple':
        return P + text + end

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
        print(form('\nOpenning gff file\n','white','bold'))
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
                        dicoRNA[id] = [Scaffold,int(start),int(end),brin]

                    elif types == 'CDS' :
                        if id not in dicoCDS.keys():
                            dicoCDS[id] = []
                        dicoCDS[id].append([Scaffold,int(start),int(end)])
    else :
        print(form("\nYou choose to don't use gff file\n",'white','bold'))
    nb = 0

    print(form('\nOpenning vcf file\n','white','bold'))

    dicoSeq = {}
    dicoInfo = {}
    listePosition = []
    test2 = True
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.CHROM not in dicoSeq.keys() :
            if int(record.POS) != 1 :
                dicoSeq[record.CHROM] = 'N'*(record.POS - 1)

        elif old_chrom == record.CHROM and old_position +1 < int(record.POS) :
            dicoSeq[record.CHROM] += 'N'*(int(record.POS) - (old_position +1 ))

        if record.INFO == {}:
            if record.CHROM not in dicoSeq.keys() :
                dicoSeq[record.CHROM] = ''
            nb+= 1


        else :
            if record.CHROM not in dicoSeq.keys():
                dicoSeq[record.CHROM] = 'N' * nb
            else :
                dicoSeq[record.CHROM] += 'N'*nb
            nb = 0
            if record.POS in listePosition:
                if len(record.REF) > len(ac):
                    # ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(record.REF)
                    ac =record.REF
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

    print(form('\nWrite fasta\n','white','bold'))
    if gff !='None' :
        with open('{0}_cds.fasta'.format(prefix), 'w') as output_cds :
            for id in sorted(dicoRNA.keys() ,key = sort_human):
                Scaffold, start, end ,brin= dicoRNA[id]
                seq = dicoSeq[Scaffold][(start-1):end]
                listeCDS = sorted(dicoCDS[id])
                cds = ''
                for elt in listeCDS:
                    cds = cds + seq[elt[1] - (int(start)):(elt[2] +1 - int(start))]
                if brin == '-':
                    seqcds = Seq(cds).reverse_complement()
                elif brin == '+':
                    seqcds = Seq(cds)
                descrip = '| position = {0}, length = {1}, brin = "{2}"'.format('{0}:{1}-{2}'.format(Scaffold,start,end), len(str(cds)), brin)
                record = SeqRecord(seqcds, id=id, name=id, description=descrip)
                SeqIO.write(record, output_cds, "fasta")

    print('\nWrite fasta Scaffold\n')
    if gff == 'None':
        with open('{0}_scaffold.fasta'.format(prefix),'w') as output_gene :

            for id in sorted(dicoSeq.keys() ,key = sort_human):
                seq = dicoSeq[id]
                record = SeqRecord(Seq(seq), id=id, name=id, description= f'| length : {len(seq)}')
                SeqIO.write(record, output_gene, "fasta")


