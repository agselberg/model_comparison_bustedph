import os 
import glob
from Bio import SeqIO

CWD = os.getcwd()
#FILE_NAMES = [os.path.basename(x).split('.fasta')[0] for x in glob.glob(os.path.join(CWD, "*.fasta"))]
##print (FA_FILES)
#for i in range(len(FILE_NAMES)):
#	lab_tree = ""
#	with open(os.path.join(CWD, "archives", FILE_NAMES[i] + ".labeled.nwk"), 'r') as file:
#		lab_tree = file.read()
#	#print(lab_tree)
#
#	species_names = []
#	for record in SeqIO.parse(os.path.join(CWD, FILE_NAMES[i] + ".fasta"), "fasta"):
#		# Extract the species name from the record's description
#		#print(record.description.split(" "))
#		#exit()
#		dummy_species= record.description.split(" ")[0]
#		species_names.append(dummy_species)
#	#print(species_names, len(species_names))
#	trait_list = []
#	#print(species_names)
#	for j in range(len(species_names)):
#		if species_names[j] not in lab_tree:
#			print("#### PROBLEM!! #### \n\n", species_names[j], "\n\n NOT IN LABELED TREE FILE \n\n", "FILE_NAMES[I]")
#			exit()
#		elif species_names[j] + '{FOREGROUND}' in lab_tree:
#			trait_list.append(species_names[j])
#	print("FILE_NAME", FILE_NAMES[i], "\n\n LEN SPECIES_NAMES", len(species_names), '\n\n LEN TRAIT_LIST', len(trait_list))
#	if len(trait_list) > 15:
#		continue
#	else:
#		with open(os.path.join(CWD, FILE_NAMES[i] + ".txt"), 'w') as f:
#			for k in range(len(trait_list)):
#				f.write(trait_list[k] + '\n')
#	#print(species_names)


species_names = []
with open(os.path.join(CWD, 'character_data.fas'), 'r') as f:
	lines = [line.strip() for line in f.readlines()]
for i in range(len(lines)):
	if lines[i] == '1':
		species_names.append(lines[i-1].lstrip('>'))

with open(os.path.join(CWD, 'housekeeping.txt'), 'w') as g:
	for i in range(len(species_names)):
		g.write(species_names[i] +'\n')
