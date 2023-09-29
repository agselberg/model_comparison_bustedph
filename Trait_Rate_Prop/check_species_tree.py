### SCRIPT TO CHECK IF LIST OF SPECIES FROM FASTA FILE MATCH SPECIES IN NEWICK TREE ###
###### THIS SCRIPT IS UNFINISHED, I STOPPED AT PRINT(FASTA_SPECIES) AND USED TIMETREE ######

nwk_file = "/home/agselberg/BUSTED-PH/busted-ph_idk/compare_methods/datasets/Trait_Rate_Prop/Tetrapoda_order.nwk"
fasta_file = "/home/agselberg/BUSTED-PH/busted-ph_idk/compare_methods/datasets/Trait_Rate_Prop/rspb20220841_si_004.fasta"

## IMPORTS ##
import re

## LOAD DATA INTO PYTHON ##
 
with open(fasta_file, 'r') as file:
	lines = file.readlines()
	fasta_species = []
	for j in range(len(lines)):
		if lines[j].startswith(">"):
			fasta_species.append(lines[j][1:].replace("\n", ""))

print(fasta_species)

with open(nwk_file, 'r') as file:
	nwk_lines = file.read()
print(nwk_lines)
for item in fasta_species:
	print(item)
species_names = []

start_index = 0
for i in range(len(nwk_lines)):
	if nwk_lines[i] == "(":
		start_index = i + 1
	elif nwk_lines[i] == "," or nwk_lines[i] == ")" or nwk_lines[i] ==':':
		species_names.append(nwk_lines[start_index:i])
		start_index = i + 1
#print(species_names)
set_nwk_sp = set(species_names)
set_fas_sp = set(fasta_species)


#print(len(fasta_species))
#print(set_nwk_sp.intersection(set_fas_sp))
#print(len(set_nwk_sp.intersection(set_fas_sp)))
