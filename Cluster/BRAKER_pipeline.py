#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package BRAKER_pipeline.py
# @author Florian Charriat

"""
	The BRAKER_pipeline script
	==========================

	:author: CHARRIAT Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to annotate with BRAKER. This script use :\n\t
	- Toogle for align multiple RNA-seq on a genome and merge the different alignments into a file\n\t
	- Samtools for file is sort\n\t  
	- Bam2hints for edit at hint format for braker\n\t
	- exonerate for align genome with annoted genome protein\n\t
	- Braker with HINTS option which use Hint of RNAseq and exonerate\n\t
	
	Be careful, please note that this script requires that filterHints.r file must be in the same directory

	Example
	-------

	>>> BRAKER_pipeline.py -d /homedir/user/work/RNAseq -o /homedir/user/work/result -c /homedir/user/work/config.txt -r /homedir/user/work/RNAseq/assembly -p homedir/user/work/annotated_genome_protein.fasta
	
	>>> BRAKER_pipeline.py -d /homedir/user/work/RNAseq -o /homedir/user/work/result -c /homedir/user/work/config.txt -r /homedir/user/work/RNAseq/assembly p homedir/user/work/annotated_genome_protein.fasta -f file.fasta
	

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display Alignement.py version number and exit
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the RNA-seq
		- \-r <path/to/reference/directory>, --refDir<path/to/reference/directory>
						path of directory that contains all the genome assembled
		- \-c <path/to/configFile>, --configFile <path/to/>configfile>
						path of the config file for toogle
		- \-p <path/to/configFile>, --proteinPath <path/to/>annotated_genome_protein.fasta>
						Path of annotated genome protein fasta file
		- \-o <path/to/output/directory>, --outDirPath <path/to/output/directory>
						path of the output directory
	Input infos for running only one genome:
		- \-f <Name/file>, --file <Name/file>
						name file genome if the user only wants to process only one genome



"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form , verifFichier , isFasta, isFastq ,  recupId



if __name__ == "__main__":

	version = '0.1'
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to align multiple RNA-seq on a genome and merge the different alignments into a file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Alignement.py version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the RNA-seq')
	filesreq.add_argument('-r', '--refDir',type = str, required=True, dest = 'refDir', help = 'Path of directory that contains all the genome assembled')
	filesreq.add_argument('-c', '--configFile',type = str, required=True, dest = 'configFile', help = 'Path of the config file for toogle')
	filesreq.add_argument('-p', '--protein',type = str, required=True, dest = 'proteinPath', help = 'Path of annotated genome protein fasta file')
	filesreq.add_argument('-o', '--outDir',type = str, required=True, dest = 'outDirPath', help = 'Path of the output directory')

	files = parser.add_argument_group('Input infos for running only one genome')
	files.add_argument('-f', '--file', type=str, required=False, dest = 'file', default ='', help = 'Name file genome if the user only wants to process only one genome')
	
	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	ref = os.path.abspath(args.refDir)
	prot = os.path.abspath(args.proteinPath)
	config = os.path.abspath(args.configFile)
	outDir= os.path.abspath(args.outDirPath)
	assembly = args.file

########### Gestion directory ##############
	directory = verifDir(directory,True)
	ref = verifDir(ref,True)
	verifFichier(config)
	verifFichier(prot)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash/'
	rna = outDir+'hints/rnaseqHints/'
	protein = outDir+'hints/proteinHints/'
	dirHints = outDir+'hints/allHints/'
	braker = outDir+'Braker/'
	name_directory = [outDir, bash,outDir+'error/',outDir+'output/',braker,rna,protein,dirHints]
	for folder in name_directory: 
		createDir(folder)

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("         Welcome in Alignement (Version " + version + ")           ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

########## Main ############################
	nbScript = 0
	run = open(outDir+'run_job_braker.sh','w')
	run.close()
	for genome in os.listdir(ref) :
		if assembly in genome and isFasta(genome):
			nbScript += 1
			IDgenome = recupId(genome)
			print("\nCréation du script d'Annotation pour : " + IDgenome)
			genomeOutDir = outDir+IDgenome
			resultMapping =	genomeOutDir+'/finalResults/'
			nameFile = bash+IDgenome+'_braker.sh'
			files = open(nameFile,'w')
			
			
################ Alignement tophat + création hints ###################
			
			#Permet de géré les sortie de sge 
			files.write('#$ -o '+outDir+'output/'+IDgenome+'.out\n#$ -e '+outDir+'error/'+IDgenome+'.err\n\n')
			
			# Permet de charger puis lancer Toogle pour un alignement
			files.write('module load bioinfo/braker/1.9\nmodule load bioinfo/exonerate/2.4.7\nmodule load bioinfo/TOGGLE/0.3.6\n')
			files.write('\ndate;\necho "Lancement script";\n\n################ Alignement tophat + création hints ###################\n\n')
			files.write('rm -r '+genomeOutDir+'\n')
			files.write('toggleGenerator.pl -d '+directory+' -r '+ref+genome+' -c '+config+' -o '+genomeOutDir+';\n')
			# Permet de récupérer tous les hits accepté pour ensuite les merger
			resultMapping = '/homedir/charriat/work/Annotation/1_tmp/align/bamMbio/scriptRes/%s/finalResults/'%(recupId(genome)) # utilisé car les alignement ont été lancé avec le script TOGGLE_ARRAY.py
			files.write('cd '+resultMapping+';\nls *.accepted_hits.bam > bamList;\n')
			
			# Permet de merger les différents hits récupérés précédement
			mergefile = 'merged_'+IDgenome+'.accepted_hits.bam'
			files.write('samtools merge -f -b bamList -c '+mergefile+';\ndate;\necho "samtools merge done";\n\n')
			
			# Permet de triée les données du fichier bam contenant tous les mapping
			sortfile  = 'merged_'+IDgenome+'.accepted_hits_sort.bam'
			files.write('java -jar /usr/local/bioinfo/picard-tools/2.7.0//picard.jar SortSam I='+mergefile+' O='+sortfile+' SORT_ORDER=coordinate;\ndate;\necho "picard-tools SortSam done";\n\n')
			
			# Permet d'utiliser l'outils bam2hints pour formater les données pour l'annotation avec Augustus ou braker
			files.write('bam2hints --minintronlen=10 --maxintronlen=1000 --maxgaplen=9 --source=M --in='+sortfile+' --out=hints_'+IDgenome+'.raw.bam;\n\ndate;\necho "bam2hints done";\n')
		
			# Permet de selectionner seulement un set de read minimum requis pour un intron avec un script R
			files.write(sys.path[0]+'/filterHints.r -s '+IDgenome+' -p '+resultMapping+'\n')
			files.write('ln -s %shints_%s.filtered.gff %s;\n'%(resultMapping,IDgenome,rna))
			


##################### Utilisation exonerate #################################

			nameOut = protein+"exonarate_"+IDgenome
			files.write('\n######################################## Exonerate #######################################\n\n')
			files.write('# Run Exonerate\nexonerate --model protein2genome --percent 95 --showtargetgff T %s %s > %s.gff3;\necho "exonerate done";\n\n'%(prot,ref+genome,nameOut))
			files.write('# Create hints\nexonerate2hints.pl --source=M --minintronlen=10 --maxintronlen=1000 --in=%s.gff3 --out=%s.hints.gff3;\ndate;\necho "exonerate2hints done";\n\n'%(nameOut,nameOut))
			
			
			
####################  Run BRAKER #######################################
			files.write('\n######################################## BRAKER #######################################\n\n')
			rnaHints = '%shints_%s.filtered.gff'%(rna,IDgenome)
			proteinHints ='%s.hints.gff3'%(nameOut)
			files.write('cat %s %s  > %sRNAseq_protein.hints_%s.gff\n'%(rnaHints,proteinHints,dirHints,IDgenome))
			files.write("awk '/intron/' %sRNAseq_protein.hints_%s.gff > %sRNAseq_protein.hints.intron_%s.gff;\n"%(dirHints,IDgenome,dirHints,IDgenome))
			hints= "%sRNAseq_protein.hints.intron_%s.gff"%(dirHints,IDgenome)
			os.system('mkdir '+braker+IDgenome) 
			files.write('\ndate;\necho "Launch Braker";\nbraker.pl --cores 24 --fungus --gff3 --species=magnaporthe_oryzae --useexisting --genome=%s --hints=%s --overwrite --alternatives-from-evidence=false --workingdir=%s\n\ndate;\necho "Job finish"'%(ref+genome,hints,braker+IDgenome))
			
 
			files.close()
			os.system('chmod 755 '+nameFile)
			run = open(outDir+'run_job_braker.sh','a')
			run.write("qsub -l mem_free=30G -N "+IDgenome+"_braker -V -q long.q '"+nameFile+"'\n")
			run.close()



############## summary message #######################
	print(form('\n-------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\tInput :')
	print('\t\t- Donnée RNAseq : '+directory[:-1])
	print('\t\t- Repertoire génome : '+ref[:-1])
	if assembly != '' :
		print('\t\t- Génome : '+assembly)
	print('\t\t- Fichier config : '+config)
	
	print('\n\tOutput :')
	print('\t\t- script bash : '+bash[:-1])
	print('\t\t- fichier a lancer : '+outDir+'run_job_braker.sh')	
	print('\t\t- Resultat des Jobs : '+ outDir[:-1])	
	
	print('\nSi vous souhaité lancer tous les '+str(nbScript)+' jobs de mapping veuillez taper la commande : ')
	print(form('\n\t\t\t\tbash '+outDir+'run_job_braker.sh\n','green','bold'))
	print(form('-------------------------------------------------------------------------------------------------------------------','red','bold'))
	
		

############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))














