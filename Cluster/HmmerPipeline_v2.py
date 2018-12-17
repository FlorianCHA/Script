#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package HmmerPiepline.py
# @author Florian Charriat

"""
	The HmmerPiepline script
	========================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 27/06/2018
	:version: 0.1

	Script description
	------------------

	This program is used to search a HMM profil for search EffecteurMAX. This program used HMMbuild et HMMsearch and mafft.

	Example
	-------

	>>> HmmerPiepline.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display HmmerPiepline.py version number and exit

	Input mandatory infos for running:
		- \-a <path/to/alignement/file>, --alignementFile <path/to/alignement/file>
						path of alignement file contrain by structure
		- \-b <path/to/sequence/database/file>, --seqdb <path/to/sequence/database/file>
						Path of sequence database
        - \-s <path/to/alignement/structure/file>, --structure <path/to/alignement/structure/file>
						Path of structure alignement file (TMalign output)

		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""

########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Import module_Flo
from module_Flo import verifDir, createDir, form, isFasta, recupId, verifFichier, fasta2dict, sort_human,indexEgale


############# Fonction ####################

def filtreSequence(file):
    """
    Remove only the redundante sequence
    """
    dico = fasta2dict(file)
    listeSequence = []
    nb = 0
    with open('%s_filtred_no_redundant.fasta'%file.replace('_filtred.fasta',''),'w') as f :
        for elt in sorted(dico.keys(),key = sort_human) :
            seq = str(dico[elt].seq)
            if seq not in listeSequence :
                nb += 1
                listeSequence.append(seq)
                sequence = Seq(seq)
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(seq)))
                SeqIO.write(record, f, "fasta")
    return nb

def filtrAln(file):
    """
    Remove sequence with no alingmement of Cystein
    """
    listeGene = []
    dico = fasta2dict(file)
    seqRef = str(dico['co39C/4-59'].seq)
    positionC = indexEgale(seqRef,'C')
    nb = 0
    with open(file,'r') as f :
        for line in f :
            if line[0] == '>':
                listeGene.append(line.split()[0].replace('>',''))

    with open(file,'w') as f ,open(file.replace('_filtred.aln','_bad.aln'),'w') as f_bad:
        for elt in listeGene :
            seq = str(dico[elt].seq)
            if seq[positionC[0]] == 'C' and seq[positionC[1]] == 'C' :
                sequence = Seq(seq)
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(seq)))
                SeqIO.write(record, f, "fasta")
            else :
                sequence = Seq(seq)
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(seq)))
                SeqIO.write(record, f_bad, "fasta")
                nb += 1
    return nb


def filtreHit(files, dico_DB):
    """
	"""
    dico = {}
    nbG = 0
    nb = 0
    with open(files, 'r') as f:
        for line in f:
            if line[0:3] == 'Mo_':
                id = line.split('/')[0]
                if id not in dico.keys():
                    dico[id] = line.split()[-1].replace('\n', '')
        # else :
        # 	dico[id] = dico[id] + line.split()[-1].replace('\n','')

    with open('%s_filtred.fasta' % files, 'w') as f_good, open('%s_remove.fasta' % files, 'w') as f_bad:

        for elt in sorted(dico.keys(), key=sort_human):
            goodAlignement = False
            index = 0
            seq = str(dico_DB[elt].seq)
            for aa in seq:
                if aa == 'C':
                    zone = seq[index + 34:index + 49]
                    if 'C' in zone:
                        sequence = Seq(seq)
                        record = SeqRecord(sequence, id=str(elt), name=str(elt), description='')
                        SeqIO.write(record, f_good, "fasta")
                        nbG += 1
                        goodAlignement = True
                        break
                index += 1
            if goodAlignement == False:
                nb += 1
                sequence = Seq(seq)
                record = SeqRecord(sequence, id=str(elt), name=str(elt), description='length : ' + str(len(sequence)))
                SeqIO.write(record, f_bad, "fasta")
    return nbG, nb # ,nbU


if __name__ == "__main__":

    version = "0.1"

    ############ Argparse #####################
    parser = argparse.ArgumentParser(prog=__file__,
                                     description='''This program is used to search a HMM profil for search EffecteurMAX. This program used HMMbuild et HMMsearch.''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help= \
        'display HmmerPiepline version number and exit')

    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-a', '--alignementFile', type=str, required=True, dest='file',
                          help='Path of alignement file contrain by structure')
    filesreq.add_argument('-s', '--structure', type=str, required=True, dest='structure',
                          help='Path of structure alignement file (TMalign output)')
    filesreq.add_argument('-b', '--seqbd', type=str, required=True, dest='db', help='Path of sequence database')
    filesreq.add_argument('-o', '--outdir', type=str, required=True, dest='outdirPath',
                          help='Path of the output directory')
    ######### Recuperation arguments ###########
    args = parser.parse_args()
    alignement = os.path.abspath(args.file)
    structure = os.path.abspath(args.structure)
    db = os.path.abspath(args.db)
    outDir = os.path.abspath(args.outdirPath)

    ########### Gestion directory ##############
    outDir = verifDir(outDir)
    name_directory = [outDir]
    createDir(name_directory)
    verifFichier(structure)

    ############### start message ########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print("\t" + form("|", 'yellow', 'bold') + form(
        "        Welcome in HmmerPiepline (Version " + version + ")         ", type='bold') + form("|", 'yellow',
                                                                                                   'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold') + '\n')

    ############## main ####################################
    profilZero = '%sprofil_0' % outDir
    os.system('hmmbuild --amino %s %s > trash ' % (profilZero, alignement))
    AlignementZero = '%salignement_0' % outDir
    newAlignement = AlignementZero
    print(form('Alignement n° 0\n', 'green', 'bold'))
    print(form('\t - Création du pattern 0 à partir de %s_filtred' % (newAlignement.split('/')[-1]), 'white', 'bold'))
    dico_fasta = fasta2dict(db)
    os.system('hmmsearch -A %s --max --nonull2 -E 1e-3 %s %s > %s' % (newAlignement, profilZero, db, '%sresult_%s' % (outDir,'0')))
    nbG, nb = filtreHit(newAlignement, dico_fasta)
    nbR = filtreSequence('%s_filtred.fasta' % newAlignement)
    # os.system('clusterProt.py -f %s_filtred.fasta -i 1 -o %s_filtred_clust.fasta --quiet' % (newAlignement, newAlignement))
    os.system('mafft --add %s_filtred_no_redundant.fasta --quiet --reorder %s > %s_filtred.aln' % (newAlignement, structure, newAlignement))
    nbC = filtrAln('%s_filtred.aln'%newAlignement)
    print(form('\t - Alignement n° 0 filtré', 'white', 'bold'))
    print(form('\t - %s : %s hits récupérés, %s hits retirés\n' % (newAlignement.split('/')[-1], nbG, nb+nbC), 'white',
               'bold'))

    i = 0
    while True:
        i += 1
        print(form('Alignement n° %s\n' % i, 'green', 'bold'))
        print(form('\t - Création du pattern %s a partir des %s hits no redondant du fichier %s_filtred.aln' % (
        i, nbR, newAlignement.split('/')[-1]), 'white', 'bold'))
        oldNb = nbG
        profil = '%sprofil_%s' % (outDir, str(i))
        # hit = '%shit_%s.fasta' % (outDir, str(i))
        os.system('hmmbuild --amino %s %s_filtred.aln > trash' % (profil, newAlignement))
        f = open(newAlignement + '_filtred.aln', 'r')
        linesOldAlignement = f.readlines()
        f.close()
        newAlignement = '%salignement_%s' % (outDir, str(i))
        print(form('\t - Utilisation du profil n° %s sur le fichier %s' % (i, db.split('/')[-1]), 'white', 'bold'))
        os.system('hmmsearch -A %s --max --nonull2 -E 1e-3 %s %s > %s' % (newAlignement, profil, db, '%sresult_%s' % (outDir, str(i))))
        nbG, nb = filtreHit(newAlignement, dico_fasta) # Filtre les hits et va ensuite crée un fasta avec toutes les séquences obtenues
        # os.system('clusterProt.py -f %s_filtred.fasta -i 1 -o %s_filtred_clust.fasta --quiet'%(newAlignement,newAlignement))
        nbR = filtreSequence('%s_filtred.fasta'%newAlignement)
        os.system('mafft --add %s_filtred_no_redundant.fasta --quiet --reorder %s > %s_filtred.aln' % (newAlignement, structure, newAlignement))
        nbC = filtrAln('%s_filtred.aln' % newAlignement)
        print(form('\t - Alignement n° %s filtré' % i, 'white', 'bold'))
        f = open(newAlignement + '_filtred.aln', 'r')
        linesNewAlignement = f.readlines()
        f.close()
        print(form('\t - %s : %s hits récupérés, %s hits retirés\n' % (newAlignement.split('/')[-1], nbG, nb+nbC), 'white',
                   'bold'))
        if linesOldAlignement == linesNewAlignement or oldNb > nbG:
            print('Il y a eu %s itération' % i)
            break

    ############## end message ###########################

    print(form("\n\t---------------------------------------------------------", 'yellow', 'bold'))
    print("\t" + form("|", 'yellow', 'bold') + form("                    End of execution                   ",
                                                    type='bold') + form("|", 'yellow', 'bold'))
    print(form("\t---------------------------------------------------------", 'yellow', 'bold'))




















