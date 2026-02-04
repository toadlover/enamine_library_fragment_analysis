#a more practical script that will run through an entire superchunk of the 2.6B enamine library and determine all fragments with counts and example ligands

#imports
import os,sys
#a basic script written with help from chatgpt to test a basic operation for pulling fragments from a single test sdf file in the 2.6B ligand library
from rdkit import Chem
from collections import defaultdict

#user input on which superchunk to investigate (from 0-530)
superchunk = int(sys.argv[1])

ROT_BOND_SMARTS = Chem.MolFromSmarts(
    "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"
)

ESTER_SMARTS = Chem.MolFromSmarts("C(=O)O[#6]")
NITRILE_SMARTS = Chem.MolFromSmarts("C#N")
#allow F and Cl
HEAVY_HALOGEN_SMARTS = Chem.MolFromSmarts("[Br,I]")

MAX_EXAMPLES = 3
MIN_HEAVY_ATOMS = 3

fragment_dict = {}

#clobber and then create a folder by the superchunk name and move into it
os.system("rm -drf " + str(superchunk))
os.system("mkdir " + str(superchunk))

os.chdir(str(superchunk))

#iterate over all chunks and subchunks within the superchunk in the library to get the ligand fragments
#iterate chunks (0-99)
for i in range(0,100):

	#derive chunk name, which is the superchunk followed by the chunk, with zeroes appended to the front to make a 5 digit string
	#i.e. superchunk 1 chunk 43 would be '00' '1' ''43' where 1 + 43 are concatenated, and then preceeding zeroes are appended
	chunk_id = str(superchunk) + str(i)

	while len(chunk_id) < 5:
		chunk_id = "0" + chunk_id

	#adding a temporary test throttle so  Ican quickly check a a few chunks quickly
	if chunk_id == "00008":
		break

	#break for at end of library
	if chunk_id == "53085":
		print("At library end and expected to look at chunk " + chunk_id)
		break

	#iterate subchunks, 0-9
	for j in range(0,10):
		#library end behavior to cut off when we hit the end and avoid moving out of bounds
		#last entry in library is superchunk 530, chunk 53084, subchunk 2
		if chunk_id == "53084" and j == 3:
		    print("At library end and expected to look at chunk " + chunk_id + " and subchunk " + str(j))
		    break

		#otherwise, copy down and extract the compressed sdf file
		print(str(superchunk) + ":" + chunk_id + ":" + str(j))

		os.system("cp /pi/summer.thyme-umw/enamine-REAL-2.6billion/" + str(superchunk) + "/" + chunk_id + "/split_new_named_" + str(j) + ".sdf.tar.gz .")
		os.system("tar -xzf split_new_named_" + str(j) + ".sdf.tar.gz")


		supplier = Chem.SDMolSupplier("split_new_named_" + str(j) + ".sdf", removeHs=True)

		for mol_idx, mol in enumerate(supplier):
		    if mol is None:
		        continue

		    #derive the ligand name as an example
		    ligand_name = (
		        mol.GetProp("_Name")
		        if mol.HasProp("_Name")
		        else f"mol_{mol_idx}"
		    )

		    #adjust the ligand name to also include the superchunk, chunk, and subchunk that it came from
		    ligand_name = str(superchunk) + ":" + chunk_id + ":" + str(j) + ":" + ligand_name


		    # Keep largest component (salt stripping)
		    frags = Chem.GetMolFrags(mol, asMols=True)
		    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

		    Chem.SanitizeMol(mol)
		    Chem.Kekulize(mol, clearAromaticFlags=False)

		    # Find rotatable bonds
		    rot_bonds = mol.GetSubstructMatches(ROT_BOND_SMARTS)

		    em = Chem.EditableMol(mol)
		    for a1, a2 in rot_bonds:
		        em.RemoveBond(a1, a2)

		    cut_mol = em.GetMol()
		    rigid_frags = Chem.GetMolFrags(
		        cut_mol, asMols=True, sanitizeFrags=True
		    )

		    seen_this_ligand = set()

		    for frag in rigid_frags:
		        # Size filter
		        if frag.GetNumHeavyAtoms() < MIN_HEAVY_ATOMS:
		            continue

		        # Chemical exclusion filters
		        if frag.HasSubstructMatch(ESTER_SMARTS):
		            continue
		        if frag.HasSubstructMatch(NITRILE_SMARTS):
		            continue
		        if frag.HasSubstructMatch(HEAVY_HALOGEN_SMARTS):
		            continue

		        frag_smiles = Chem.MolToSmiles(
		            frag,
		            canonical=True,
		            isomericSmiles=False
		        )

		        if frag_smiles in seen_this_ligand:
		            continue
		        seen_this_ligand.add(frag_smiles)

		        if frag_smiles not in fragment_dict:
		            fragment_dict[frag_smiles] = {
		                "count": 1,
		                "examples": [ligand_name]
		            }
		        else:
		            fragment_dict[frag_smiles]["count"] += 1
		            if len(fragment_dict[frag_smiles]["examples"]) < MAX_EXAMPLES:
		                fragment_dict[frag_smiles]["examples"].append(ligand_name)

		os.system("rm split_new_named_" + str(j) + ".sdf.tar.gz")
		os.system("rm split_new_named_" + str(j) + ".sdf")


#output the fragment dictionary to a csv file
write_file = open(str(superchunk) + "_fragments.csv", "w")
write_file.write("fragment_smiles,occurences,example_ligands\n")
for frag in fragment_dict.keys():
	write_file.write(frag + "," + str(fragment_dict[frag]["count"]) + ",")
	#write the examples
	for ligname in fragment_dict[frag]["examples"]:
		write_file.write(str(ligname) + ";")

	write_file.write("\n")
