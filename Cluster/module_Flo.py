#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package module_Flo.py
# @author Florian Charriat
#__docformat__ = "restructuredtext en"

"""
	The module_Flo module
	=====================

	author: CHARRIAT Florian\n
	contact: florian.charriat@inra.fr\n
	date: 21/03/2018\n
	version: 0.1\n

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

def createDir(directory):
	''' Permet de vérifier si un dossier existe, si ce n'est pas le cas, 
	le dossier sera crée
	'''
	if not os.path.exists(directory):
		 	os.makedirs(directory)
	return


def verifFichier(fichier):
	''' Permet de vérifier si un fichier existe
	
	'''
	if os.path.exists(fichier):
		return
	else :
		raise ValueError("ERROR the file '%s' doesn't exist, please check if your files exists" % fichier)	







def verifDir(directory,check = False):
	'''
	Permet de mettre en forme le chemin du dossier pour être utilisé dans un script, 
	la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute
	'''
	if directory.endswith('/') == False :
		directory = directory +'/'
	if check :
		if os.path.isdir(directory):
			return directory
		else :
			raise ValueError("ERROR the directory'%s' is not valid path, please check if your directory exists" % directory)
	else :
		return directory

def form(text,col = 'white' ,type = 'none') :
	'''
	Permet de mettre en forme les textes afficher sur le terminale,
	text : La chaine de caractère a transformer
	col : La couleur souhaité entre les couleurs red, green,yellow,orange,blue et purple
	type : formatage du texte souhaité entre les formats bold, underline, blind et highligth.
	 Si plusieurs type sont souhaité il est possible de faire une liste (exemple, type = ['bold','underline'])
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











