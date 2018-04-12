#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package module_Flo.py
# @author Florian Charriat
#__docformat__ = "restructuredtext en"

"""
	The module_Flo module
	=====================

	:author: CHARRIAT Florian\n
	:contact: florian.charriat@inra.fr\n
	:date: 21/03/2018\n
	:version: 0.1\n

	Use it to import very handy functions.

	Example:

	>>> from module_Flor import createDir
	>>> createDir('resultat')
	
"""
##################################################
## Modules
##################################################
## Python modules
import argparse, os



####### FUNCTION ################

################################ Fonction repertoire ##################################################"

def createDir(Listedirectory):
	'''Permet de vérifier si un dossier existe, si ce n'est pas le cas, 
	le dossier sera crée.
	
	:Parameters:
	     Listedirectory
		liste de repertoire à créer
	'''
	
	if type(Listedirectory) != list:
		if not os.path.exists(Listedirectory):
			os.makedirs(Listedirectory)
	
	else :
		for directory in Listedirectory:
			if not os.path.exists(directory):
				 	os.makedirs(directory)
	return
	


def verifDir(directory,check = False):
	'''
	Permet de mettre en forme le chemin du dossier pour être utilisé dans un script,la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute. La fonction peut aussi verifier qu'un repertoire existe.
	
	:Parameters:
	     directory
		Path du dossier
	     check : bool
	     	Si check = True, la fonction va aussi verifier que le repertoire existe
	     	
	'''
	if directory.endswith('/') == False :
		directory = directory +'/'
	if check :
		if os.path.isdir(directory):
			return directory
		else :
			raise ValueError(form("ERROR the directory '%s' is not valid path, please check if your directory exists" % directory,'red','bold'))
	else :
		return directory
		
################################## Fonction fichier ############################################"

def verifFichier(fichier):
	'''Permet de vérifier si un fichier existe.
	
	:Parameters:
	     fichier
		Path du fichier
	    
	'''
	if os.path.exists(fichier):
		return
	else :
		raise ValueError(form("ERROR the file '%s' doesn't exist, please check if your files exists" % fichier,'red','bold'))
		
##################################### Fonction fichier fasta/fastq #################################################"

def isFasta(fichier):
	'''Permet de vérifier si un fichier est au format fasta	
	
	:Parameters:
	     fichier
		Path du fichier	
	'''
	if fichier.endswith('.fasta') or fichier.endswith('.fa') or fichier.endswith('.fasta.gz') or fichier.endswith('.fa.gz'):
		return True
	else :
		return False
		

def isFastq(fichier):
	'''Permet de vérifier si un fichier est au format fastq	
	
	:Parameters:
	     fichier
		Path du fichier
	'''
	if fichier.endswith('.fastq') or fichier.endswith('.fq') or fichier.endswith('.fastq.gz') or fichier.endswith('.fq.gz'):
		return True
	else :
		return False
		
def recupId(fichier):
	'''Permet de récupéré le nom du fichier sans l'extension fasta ou fastq
	
	:Parameters:
	     fichier
		Path du fichier
	'''
	# Traitement pour fichier fasta
	fichier = fichier.replace('.fasta.gz','')
	fichier = fichier.replace('.fa.gz','')
	fichier = fichier.replace('.fasta','')
	fichier = fichier.replace('.fa','')
	
	# Traitement pour fichier fastq
	fichier = fichier.replace('.fastq.gz','')
	fichier = fichier.replace('.fq.gz','')
	fichier = fichier.replace('.fastq','')
	fichier = fichier.replace('.fq','')
	
	return fichier 
	
	
	
#################################### Fontion formatage texte ################################################"

def form(text,col = 'white' ,type = 'none') :
	'''
	Permet de mettre en forme les textes afficher sur le terminale.
		
	:Parameters:
	     text: string
		Le texte à transformer
	     col: str
		La couleur souhaité entre les couleurs red, green,yellow,orange,blue et purple
	     text: str
		 str ou liste de str du format à appliquer (bold, underline, blind et highligth)
	'''
	W  = '\033[0'  # white (normal)
	R  = '\033[31' # red
	G  = '\033[32' # green
	Y  = '\033[33' # yellow
	O  = '\033[33' # orange
	B  = '\033[34' # blue
	P  = '\033[35' # purple
	end = '\033[0m'  # white (normal)
	Bold = ';1'
	underline = ';4'
	blind = ';5'
	highlight =';7' 	
	text = 'm'+text		
	if 'bold' in type :
		text = Bold + text
	if 'underline' in type :
		text = underline + text
	if 'highlight' in type :
		text = blind + text
	if 'highlight' in type :
		text = highlight + text
	if col == 'red' :
		return R+text+end
	elif col == 'white' :
		return W+text+end
	elif col == 'green' :
		return G+text+end
	elif col == 'yellow' :
		return Y+text+end
	elif col == 'orange' :
		return O+text+end
	elif col == 'blue' :
		return B+text+end
	elif col == 'purple' : 
		return P+text+end











