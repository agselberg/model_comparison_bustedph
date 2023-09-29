########################################################
# converting coevol data from 2012 paper to hyphy format
# by Avery Selberg
########################################################

from Bio import SeqIO
import os

cwd = os.getcwd()

trait_list = '/home/agselberg/BUSTED-PH/busted-ph_practice/compare_methods/datasets/coevol/data/plac/plac.lht'
fasta_plac = '/home/agselberg/BUSTED-PH/busted-ph_practice/compare_methods/datasets/COEVOL_2012/fastas/ZFX_73placentalsLongevity.fasta'
fasta_mam = '/home/agselberg/BUSTED-PH/busted-ph_practice/compare_methods/datasets/COEVOL_2012/fastas/ZFX_78mammalsLongevity.fasta'



### LOAD FASTA DATA INTO PYTHON ###
plac_species = set()
mam_species = set()

# iterate over the sequences in each file and extract the species name from the header
for record in SeqIO.parse(fasta_plac, "fasta"):
    species_name = record.description # assumes species name is the second word in the header
    plac_species.add(species_name)

for record in SeqIO.parse(fasta_mam, "fasta"):
    species_name = record.description # assumes species name is the second word in the header
    mam_species.add(species_name)
plac_species = list(plac_species)
mam_species = list(mam_species)


#fasta_species = {'plac_species' : plac_species, 'mam_species' : mam_species}
fasta_species = {'plac_species' : plac_species}
#print(fasta_files)


# print the sets of species names
#print("Species in file1: ", plac_species)
#print("Species in file2: ", mam_species)

#exit()

### LOAD TRAIT DATA INTO PYTHON ###
 ## TRAIT DATA CONSISTS OF A MATRIX OF FOUR COLUMNS ##


## LOAD IN ALL DATA AS A LIST
traits = []
with open(trait_list, 'r') as f:
    next(f)
    headers = next(f).strip().split()
    for line in f:
        traits.append(line)


## CONVERT LIST OF TRAITS, HEADERS TO 'LIST OF LISTS' AND CHANGE DATA TYPES ##
headers = headers[2:]

traits = [t.strip().split('\t') for t in traits]


new_traits = []
for trait in traits:
    new_trait = []
    for item in trait:
        if '.' in item or (item.startswith('-') and item[1:].isdigit()):
            new_trait.append(float(item))
        elif item.isdigit():
            new_trait.append(int(item))
        else:
            new_trait.append(item)
    new_traits.append(new_trait)
traits = new_traits


#print(traits)
#exit()

### DOUBLE CHECK ALL SPECIES HAVE A TRAIT, FOR EACH TRAIT ###
traits_tuple = [(trait[0], [trait[1], trait[2], trait[3]]) for trait in traits]
traits_dict = {}
#print(traits_tuple)
#exit()
for j in range(len(headers)):
    traits_dict[headers[j]] = sorted(traits_tuple, key=lambda x: x[1][j])
    #print(traits_dict)
    species_list_traits = [key[0] for key in traits_dict[headers[j]]]
    #print('SPECIES_ LIST', species_list)
    #exit()
    species_no_traits = []
    species_trait_sp = []
    for k in fasta_species:
        print(fasta_species[k])
        for i in range(len(fasta_species[k])):
               if fasta_species[k][i] not in species_list_traits:
                       species_no_traits.append(fasta_species[k][i])
               elif fasta_species[k][i] in species_list_traits:
                       species_trait_sp.append(fasta_species[k][i])

        print("SPECIES NO TRAITS", headers[j])
        print(species_no_traits)
        if len(species_no_traits) > 0:
            print('SPECIES WITHOUT TRAITS! EXITING NOW')
            print('SPECIES', k, 'TRAITS', headers[j])
            exit()
        print("SPECIES WITH TRAITS", headers[j])
        print(species_trait_sp)
        #exit()
        ### TRIM TRAITS TO TOP 20% ##
        start_index = int(len(traits_dict[headers[j]]) * 0.8)
        traits_dict[headers[j]] = traits_dict[headers[j]][start_index:]
        #print('TRAITS_DICT', traits_dict)
        species_to_write = [key[0] for key in traits_dict[headers[j]]]
        #print('TRAIT NAME', headers[j])
        #print('SPECIES_TO_WRITE', species_to_write)
        #exit()
        print('OUTPUT LOCATION', os.path.join(cwd, k + headers[j] + '.txt'))
        with open(os.path.join(cwd, k + headers[j] + '.txt'), 'w') as txt_file:
            for species in species_to_write:
                txt_file.write(species + '\n')

