


f = open('/homedir/gladieux/work/magMax_project/6_EffecteurMax/2_hmmer/v1/test_pipeline_on_Jerone_alignement/alignement_13_filtred','r')
lines = f.readlines()
f.close()
effecteurMax = []
for line in lines :
	if line[0] != '#' and line[0] != '\n' :
		effecteurMax.append(line.split('/')[0])

	


#f = open('/homedir/gladieux/work/magMax_project/4_Orthologie/0_rawdata/Results_Jul10_1/Orthogroups.GeneCount.csv','r')
#lines = f.readlines()
#f.close()
#core = open('/homedir/gladieux/work/magMax_project/4_Orthologie/coreGenome.txt','w')
#core.write(lines[0])
#entete = lines[0].split()
#lines = lines[1:]
#dico_count = {}
#for line in lines :
#	groupe = line.split()[0]
#	count = line.split()[1:-1]
#	dico_count[groupe]=count
	#a = groupe
	#rice = 0
	#other = 0
	#liste_other = []
	#if '0' not in count :
		#core.write(line)
	
		
#core.close()	
	
f = open('/homedir/gladieux/work/magMax_project/4_Orthologie/0_rawdata/Results_Jul10_1/Orthogroups.csv','r')
lines = f.readlines()
f.close()
lines = lines[1:]
for line in lines :
	groupe = line.split()[0]
	count = line.split()[1:]
	a = groupe
	for elt in count :
		if elt in effecteurMax :
			print(groupe, len(count))
			print(elt)
			break

		
		
		
		
		
		
		
		














def indexEgale(liste,str):
	index = []
	for i,e in enumerate(liste):
		if e == str :
			index.append(i)
	return index
	
def indexDif(liste,str):
	index = []
	for i,e in enumerate(liste):
		if e != str :
			index.append(i)
	return index

