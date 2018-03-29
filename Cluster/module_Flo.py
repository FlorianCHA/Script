#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package module_Flo.py
# @author Florian Charriat
#__docformat__ = "restructuredtext en"
"""
	The module_Flo module
	======================
	author: CHARRIAT Florian
	contact: florian.charriat@inra.fr
	date: 21/03/2018
	version: 0.1

	Use it to import very handy functions.

	Example:

	>>> from MODULES_SEB import dict2txt
	>>> dico = {"key1":"value1","key2":"value2","key3":"value3"}
	>>> dict2txt(dico)
	key1	value1
	key2	value2
	key3	value3
	
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



def verifDir(directory):
	'''
	Permet de mettre en forme le chemin du dossier pour être utilisé dans un script, 
	la fonction vérifie si il y a bien un '/' à la fin du chemin, sinon il le rajoute
	'''
	if directory.endswith('/') == False :
		directory = directory +'/'

	return directory
