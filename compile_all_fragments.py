#this will go through and compile all fragments from all csv files

#imports
import os,sys
from rdkit import Chem

#create a master dictionary that holds all of the unique fragments with counts and examples
fragment_dict = {}

#iterate over each */*_fragments.csv file
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith("_fragments.csv"):
			print(r + "/" + file)

			#read the file
			read_file = open(r + "/" + file, "r")

			for line in read_file.readlines():
				#skip the header
				line_stripped = line.strip()

				if line.startswith("fragment_smiles,occurences,example_ligands"):
					continue

				#pull the line data
				smiles_str = line_stripped.split(",")[0]
				count = int(line_stripped.split(",")[1])
				#exist as list, separated by semicolons
				examples = line_stripped.split(",")[2].split(";")

				#add the fragment to the dictionary or update the existing fragment entry
				#get a rdkit moleucle and convert to a rdkit-canon smiles string for comparison (probably unecessary given the original strings came from rdkit, but safe)
				mol = Chem.MolFromSmiles(smiles_str)
				canon = Chem.MolToSmiles(mol, canonical=True)

				#see if the fragment smiles is in the dictionary
				#not in
				if canon not in fragment_dict:
					fragment_dict[canon] = {
						"count": count,
						"examples": examples
					}
				#in
				else:
					fragment_dict[canon]["count"] += count

					#add examples until we have 3 or have none to see
					if len(fragment_dict[canon]["examples"]) < 3:
						for i in range(0,len(examples)):
							if len(fragment_dict[canon]["examples"]) < 3:
								fragment_dict[canon]["examples"].append(examples[i])

#output the fragment dictionary to a csv file
write_file = open("all_fragments.csv", "w")
write_file.write("fragment_smiles,occurences,example_ligands\n")
for frag in fragment_dict.keys():
	write_file.write(frag + "," + str(fragment_dict[frag]["count"]) + ",")
	#write the examples
	for ligname in fragment_dict[frag]["examples"]:
		write_file.write(str(ligname) + ";")

	write_file.write("\n")
				