#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package scriptName.py
# @author Florian CHARRIAT

## Python modules
import argparse, os , glob, re
## BioPython
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def sort_human(s, _nsre=re.compile('([0-9]+)')):
	""" Sort the list in the way that humans expect, use list.sort(key=sort_human) or sorted(list, key=sort_human)).
	"""
	try:
		return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]
	except TypeError:
		if not isinstance(s, int):
			print("WARNNING MODULES_SEB::sort_human : List %s value not understand so don't sort \n" % s)


	return s

def fasta2dict(filename):
	"""
	Function that take a file name (fasta), and return a dictionnary of sequence
	"""
	with open(filename, "rU") as fastaFile:
		return SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))


pathFasta = 'OG.fasta'
dico = fasta2dict(pathFasta)
nb = 0
with open('OG_filter.fasta','w') as output_handle :
	for elt in sorted(dico.keys(),key = sort_human):
		seq = str(dico[elt].seq)
		if seq[0] == 'M' and  '*' not in seq[0:len(seq)-1] and 'X' not in seq and len(seq) > 20:
			if seq[-1] != '*' :
				seq = seq +'*'
			sequence = Seq(seq)
			record = SeqRecord(sequence, id=elt, name=elt, description=dico[elt].description)
			SeqIO.write(record, output_handle, "fasta")
		else :
			nb += 1

print(f'{nb} séquences supprimées')