blast_output = snakemake.input[0]
orthologueCount = snakemake.input[1]
output_tmp = snakemake.output[0]
OG = snakemake.params[0] # Option a voir
dico_OG = dico_OG(orthologueCount)
f = open(blast_output,'r')
lines = f.readlines()
f.close()
f = open('%scorrection/%s.fasta'%(output_tmp),'w')
liste = []
startFile = True
for line in lines :
	if line[0:6] == 'Query=' :
		lenQuery = line.split()[-1]
	elif line[0] == '>' and startFile == True:
		lineSplit = line.split()
		name = lineSplit[1].split('_')[0]
		Num_scaffold = lineSplit[1].split('_')[-1]
		start = True
		startFile = False
	elif line[0] == '>' :
		if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
			if name not in dico_OG[OG] :
				dico_OG[OG].append(name) # A voir
				start,end,sens = functionSens(pos1,pos2)
				sequence = Seq(seq)
				record = SeqRecord(sequence,id=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)),name=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)), description= 'length : '+ str(len(sequence)))
				SeqIO.write(record,f, "fasta")
				nbC += 1

		lineSplit = line.split()
		name = lineSplit[1].split('_')[0]
		Num_scaffold = lineSplit[1].split('_')[-1]
		start = True
	elif 'Identities' in line :
		lineSplit = line.split()
		PcIdent = int(lineSplit[3].replace('(','').replace('%),',''))
		couverture = int(lineSplit[2].split('/')[-1])/int(lenQuery) *100
	elif 'Sbjct' in line and start == True:
		lineSplit = line.split()
		pos1 = lineSplit[1]
		startAA = lineSplit[2][0]
		seq = lineSplit[2]
		start = False
	elif 'Sbjct' in line and start == False:
		lineSplit = line.split()
		pos2 = lineSplit[1]
		endAA = lineSplit[2][-1]
		seq = seq +  lineSplit[2]
	
if startAA == 'M' and endAA == '*' and couverture == 100 and PcIdent > 85 :
	if name not in dico_OG[OG] :
		start,end,sens = functionSens(pos1,pos2)
		record = SeqRecord(sequence,id=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)),name=str('%s_putative_%s:%s-%s'%(name,Num_scaffold,start,end)), description= 'length : '+ str(len(sequence)))
		SeqIO.write(record,f, "fasta")
		nbC += 1
f.close()

