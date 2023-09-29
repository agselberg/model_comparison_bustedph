
########################################################
## Script to unconcatenate coevol files from 2012 paper
## By Avery Selberg
########################################################


# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## IMPORTANT ## DON'T USE 78 MAMMALS FOR THIS ANALYSES,
## HAVEN'T OBTAINED TRAIT INFROMATION ##

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# imports
import os


# load in files, set wd

cwd = os.getcwd()
nex_files = [f for f in os.listdir(cwd) if f.endswith('.nex')]

for i in range(len(nex_files)):
	tree = nex_files[i].split('.nex')[0] + '_intree.tre'
	tree = tree.replace('&', r'\&')
	#print(tree)
	#exit()
	nex = nex_files[i].split('Evolution')[1].split('.nex')[0]
	#print(nex)
	with open(nex_files[i], 'r') as f:
		lines = f.readlines()
	#print(nex_files[i])
	filtered_lines = [s for s in lines if s != '\n']
	source_lines = [(i, line) for i, line in enumerate(filtered_lines) if line.startswith('[Source of sequence data:')]
	end_lines = [l for l, line in enumerate(filtered_lines) if line == 'end;\n'][0]
	for j in range(len(source_lines)):	
		file_lines = (0,0)
		if (j != len(source_lines) -1 ):
			file_lines = range(source_lines[j][0]+1, source_lines[j+1][0])
		if ( j == len(source_lines)-1):
			file_lines = range(source_lines[j][0] + 1, end_lines-1)
		fasta_list = []
		count = 0
		for k in file_lines:
			fasta_list.append(filtered_lines[k].split('\t'))
			fasta_list[count] = [s.strip() for s in fasta_list[count] if s.strip()] 
			count += 1
		fasta_name = source_lines[j][1].split('data:')[1].split('_gblocks')[0].strip() + nex 
		trait_names = ['Mass', 'Maturity', 'Longevity']
		for k in range(len(trait_names)):
		#	with open(os.path.join(cwd, 'fastas', fasta_name + trait_names[k] + '.fasta'), 'w') as s:
		#		for l in range(len(fasta_list)):
		#			s.write('>' + fasta_list[l][0] + '\n')
		#			s.write(fasta_list[l][1] + '\n')

		#	#os.system('cp plac_' + trait_names[k] + '_noStop.txt ' + os.path.join(cwd, 'traits', fasta_name + trait_names[k] + '.txt'))
		#	#os.system('cp ' + tree + ' ' + os.path.join(cwd, 'trees', fasta_name + trait_names[k] + '.nwk')
			os.system('cp plac_species' + trait_names[k] + '.txt ' + os.path.join(cwd, 'traits', fasta_name + trait_names[k] + '.txt'))

			

