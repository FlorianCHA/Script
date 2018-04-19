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

	>>> from module_Flo import createDir
	>>> createDir('resultat')
	
"""
##################################################
## Modules
##################################################
## Python modules
import argparse, os , glob



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
	'''Permet de vérifier si un fichier est au format fasta, renvoie True si le fichier est au format fasta.
	
	:Parameters:
	     fichier
		Path du fichier	
	'''
	if fichier.endswith('.fasta') or fichier.endswith('.fa') or fichier.endswith('.fasta.gz') or fichier.endswith('.fa.gz'):
		return True
	else :
		return False
		

def isFastq(fichier):
	'''Permet de vérifier si un fichier est au format fastq, renvoie True si le fichier est au format fastq.
	
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
	     text 
		Le texte à transformer
	     col 
		La couleur souhaité entre les couleurs red, green, yellow, orange, blue et purple
	     text
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
		
		
################ Class directory #################

class directory(str):
	'''
	Class which derives from string.
	Checks that the string is and path to valid directory and not empty
	'''
	def __init__(self, path = None):
		'''
		Initialise les variables
		'''
		self.listFiles = []
		self.listDir = []
		
		#Change le chemin en chemin absolu
		self.path = os.path.abspath(path)
		
		#Appel fonction pour avoir les informations voulu
		self.exist() # Permet de verifier si le repertoire existe
		self.verif() # Permet de verifier si le chemin est correctement ecrit
		self.listAll() # Permet de récupérer tous ce qui se trouve dans le répertoire
		self.type() # Permet de crée les deux listes, celle contenant que les fichiers et l'autre contenant que les dossiers
		
	
	def __str__(self):
		"""Fonction qui permet de formater le text de sortie lors du print"""
		return """
\033[32;1mpathDirectory\033[0m : %s\n
\033[32;1mlistPath\033[0m : %s\n
\033[32;1mlistDir\033[0m : %s\n
\033[32;1mlistFiles\033[0m : %s\n
""" % (self.path, str(self.listAll), str(self.listDir), str(self.listFiles))

	def exist(self):
		'''
		Fonction qui vérifie que le repertoire existe
		'''
		if os.path.isdir(self.path) != True :
			raise ValueError(form("ERROR the path '%s' is not valid path, please check if your directory exists" % self.path,'red','bold'))
								
	def verif(self):
		'''
		Permet de mettre en forme le chemin du dossier pour être utilisé dans un script,la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute. La fonction peut aussi verifier qu'un repertoire existe.     	
		'''
		if self.path.endswith('/') == False :
			self.path = self.path +'/'
		
	def listAll(self):
		'''
		liste tous ce qui se trouve dans le repertoire
		'''
		self.listAll = glob.glob(self.path+"*")
	def type(self):
		'''
		Créé une liste de fichier et une liste de répertoire
		'''
		for elt in self.listAll :
			if os.path.isdir(elt) == True :
				self.listDir.append(elt)
			elif os.path.exists(elt) == True :
				self.listFiles.append(elt)
			else :
				print(form('Attention, le repertoire ne contients pas que des fichiers et des sous-repertoires','red','bold'))
	
	def listExt(self,extension):
		'''
		Permet de créer une liste de fichier d'extension donnée
		:Parameters:
	   	     extension
	   	     	extension des fichiers à lister
	   	'''
		if extension == 'fasta':
			extension = ('fasta','fa')
		if extension == 'fastq':
			extension = ('fq','fastq')
		if extension == 'allfasta':
			extension = ('fasta','fa','fasta.gz','fa.gz')
		if extension == 'allfastq':
			extension = ('fq','fastq','fastq.gz','fq.gz')
		listeExt = []
		for elt in self.listAll :
			split = elt.split('/')
			if '.' in split[-1] :
				extensionElt = split[-1].split('.')[-1]
			else :
				extensionElt = 'directory'
			if extensionElt in extension :
				listeExt.append(elt)
		if len(listeExt) == 0 :
			raise ValueError(form("ERROR, the path '%s' doen't contain %s files , please check if your directory exists" % (self.path,extension),'red','bold'))
		return listeExt
			










