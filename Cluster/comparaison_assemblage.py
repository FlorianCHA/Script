#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-



########## Module ###############
## Python modules
import argparse, os, sys

#Import MODULES_SEB
from module_Flo import verifDir, createDir



if __name__ == "__main__":

	version = '0.1'
	
########### Gestion directory ##############
directory = '/gs7k1/projects/BGPI/becphy/datas/1_temporary/assemblies/AllAssemblyFasta-Seb-paired/'
output = '/homedir/charriat/work/data_seb.txt'
result = open(output,'w')
result.write('Id\tNum\tN50\tL50\n')
result.close()
for files in os.listdir(directory):
	if files.endswith('.fa.N50'):
		data = open(directory+files,'r')
		lines = data.readlines()
		result = open(output,'a')
		result.write(files.replace('.fa.N50','')+'\t'+lines[0].split('\t')[1][:-1]+'\t'+lines[3].split('\t')[1][:-1]+'\t'+lines[4].split('\t')[1])
		result.close
