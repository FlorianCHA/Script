#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package gvcfSearch.py
# @author Florian Charriat



########## Module ###############
## Python modules
import argparse, os, sys,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


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

def recupAC(line,MQmin,DPmin,QDmin):
    '''
    Function that allows to recover the nucleic acid according to the quality of mapping
    '''
    IUPAC = {"AT": "W", "CG": "S", "AC": "M", "GT": "K", "AG": "R", "CT": "Y","TA": "W", "GC": "S", "CA": "M", "TG": "K", "GA": "R", "TC": "Y"}
    K, pos, _, REF, ALT, QUAL, _, info_map, _, info = line.split()
    if ',<NON_REF>' in ALT :
        ALT = ALT.replace(',<NON_REF>','')
    ALT = ALT.split(',')
    MQ = float(info_map.split('MQ=')[1].split(';')[0])
    DP = int(info_map.split('DP=')[1].split(';')[0])
    if diploidy:
        infoDiploidy = info.split(':')[0]
    if DP < DPmin or MQ < MQmin:
        ac = 'N'*len(REF)
    elif ALT[0] == None :
        ac = REF
    elif ALT[0] != None and DP >= DPmin and MQ >= MQmin:
        if float(QUAL) >= QDmin:
            if diploidy :
                allele1 = int(infoDiploidy.split('/')[0])
                allele2 = int(infoDiploidy.split('/')[1])
                if allele1 == 0 and allele2 == 0 :
                    ac = REF
                if allele1 == 0 or allele2 == 0 :
                    if len(REF) == 1 and len(ALT[allele1-1]) == 1 and len(ALT[allele1-1]) == 1:
                        ac = IUPAC[str(REF) + str(ALT[allele1-1])]
                    else :
                        ac = 'N'*len(REF)

                elif allele1 > 0 and allele2 > 0 and allele1 == allele2 :
                    if allele1 > len(ALT) :
                        ac = 'N' * len(REF)
                        print(infoDiploidy)
                        print(ALT)
                        print(K,pos)

                    elif len(REF) == 1 and len(ALT[allele1-1]) == 1:
                        ac = ALT[allele1-1]
                    elif len(REF) >= 2 and len(ALT[allele1-1]) == 1:
                        ac = str(ALT[allele1-1]) + 'N' * (len(REF) - 1)
                    elif len(ALT[allele1-1]) >= 2:
                        ac = 'N' *len(REF)


                elif allele1 > 0 and allele2 > 0 and allele1 != allele2 :
                    if allele1 > len(ALT) or  allele2 > len(ALT):
                        ac = 'N' * len(REF)
                        print(infoDiploidy)
                        print(ALT)
                        print(K, pos)

                    elif len(REF) == 1 and len(ALT[allele1-1]) == 1 and len(ALT[allele2-1]) == 1:
                        ac = IUPAC[str(ALT[allele1-1]) + str(ALT[allele2-1])]
                    else :
                        ac = 'N'*len(REF)
            else :
                ac = ALT[0]
        else:
            ac = 'N'*len(REF)
            # print("Warning, the gene {} have a SNP (QD = {}) which don't pass the filter (QD < {})".format(id,QD,QDmin))
    return str(ac)

if __name__ == "__main__":
    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to search a list of gene in a vcf file and give in output
                                     CDS,protein and gene fasta file.''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display gvcfSearch version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-f', '--file', type=str, required=True, dest='vcf_file',
                          help='Path of vcf file which be used')
    filesreq.add_argument('-a', '--assembly', type=str, required=True, dest='fasta',
                          help='Path of assembly fasta file which be used as reference')

    filesreq.add_argument('-g', '--gff', type=str, required=False,default = 'None', dest='gff',
                          help='Path of the gff file path of the reference genome used to create the vcf')

    files = parser.add_argument_group('Input infos for running with default values')
    files.add_argument('-MQ', '--MQ', type=int, required=False, default=20, dest='MQmin',
                       help='The minimum mapping quality supporting to accept the mapping (default = 20)')

    files.add_argument('-DP', '--DP', type=int, required=False, default= 10, dest='DPmin',
                          help=' The minimum sequencing depth supporting to accept the mapping ( default = 10)')
    files.add_argument('-Q', '--QUAL', type=int, required=False, default= 1, dest='QDmin',
                          help='The minimum QUAL value to accept the variant ( default = 1)')
    files.add_argument('-diploidy', '--diploidy',action='store_true', dest='diploidy',
                          help='Use this option if you are using a diploid organism')
    files.add_argument('-p', '--prefix', type=str, required=False, default='gvcfSearch', dest='prefix',
                          help='prefix for the output file')

    ######### Recuperation arguments ###########
    args = parser.parse_args()
    vcf_file = os.path.abspath(args.vcf_file)
    gff = args.gff
    fasta = os.path.abspath(args.fasta)
    prefix = args.prefix
    MQmin = args.MQmin
    QDmin = args.QDmin
    DPmin = args.DPmin
    diploidy = args.diploidy
    ########### Gestion directory ##############
    verifFichier(vcf_file)
    verifFichier(fasta)
    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print(
        "\t" + form("|", 'yellow', 'bold') + form("            Welcome in gvcfSearch (Version " + version + ")         ",
                                                  type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

    ############### Main ########################

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

                if types == 'mRNA' :
                    id = id.split(';')[0].split('=')[1]

                    dicoRNA[id] = [Scaffold,int(start),int(end),brin]

                elif types == 'CDS' :
                    id = id.split('Parent=')[-1]
                    if id not in dicoCDS.keys():
                        dicoCDS[id] = []
                    dicoCDS[id].append([Scaffold,int(start),int(end)])
    else :
        print(form("\nYou choose to don't use gff file\n",'white','bold'))
    nb = 0
    print(form('Openning fasta file\n', 'white', 'bold'))
    dico_assembly = {}
    dico  = fasta2dict(fasta)
    for elt in dico.keys() :
        dico_assembly[elt] = str(dico[elt].seq)
    print(form('\nOpenning vcf file\n','white','bold'))
    with open(vcf_file,'r') as vcf_read :
        for line in vcf_read :
            if line[0] != '#' :
                if 'END=' in line :
                    K, pos, _, REF, ALT, _, _, end, _, info = line.split()
                    DP_min_block = int(info.split(':')[3])
                    if DP_min_block < DPmin or DP_min_block == 0 :
                        end = int(end.replace('END=',''))
                        start = int(pos)
                        if start == 1 :
                            dico_assembly[K] = 'N' * (end - start) + dico_assembly[K][end - 1:len(dico_assembly[K])]
                        else :
                            dico_assembly[K] = dico_assembly[K][0:(start-1)] +'N'*(end-start) + dico_assembly[K][end-1:len(dico_assembly[K])]
                else :
                    K, pos, _, REF, ALT, QUAL, _, info_map, _, info = line.split()
                    ac = recupAC(line, MQmin, DPmin, QDmin)
                    pos = int(pos)
                    dico_assembly[K] = dico_assembly[K][0:pos-1] +ac + dico_assembly[K][pos+(len(ac)-1):len(dico_assembly[K])]


    print(form('\nWrite fasta\n','white','bold'))
    if gff !='None' :
        with open('{0}_cds.fasta'.format(prefix), 'w') as output_cds :
            for id in sorted(dicoRNA.keys() ,key = sort_human):
                Scaffold, start, end ,brin= dicoRNA[id]
                seq = dico_assembly[Scaffold][(start-1):end]
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
            for id in sorted(dico_assembly.keys() ,key = sort_human):
                seq = dico_assembly[id]
                record = SeqRecord(Seq(seq), id=id, name=id, description= f'| length : {len(seq)}')
                SeqIO.write(record, output_gene, "fasta")


