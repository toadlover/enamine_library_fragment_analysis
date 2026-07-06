#a more practical script that will run through an entire superchunk of the 2.6B enamine library and determine all fragments with counts and example ligands

#imports
import os,sys
#a basic script written with help from chatgpt to test a basic operation for pulling fragments from a single test sdf file in the 2.6B ligand library
from rdkit import Chem
from collections import defaultdict

#user input on which superchunk to investigate (from 0-530)
superchunk = int(sys.argv[1])

#definition of the non-rotatable fragment wildcard for fragment grabbing
ROT_BOND_SMARTS = Chem.MolFromSmarts(
    "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"
)

#remove these filters for curation from the 64b library because we want everything
ESTER_SMARTS = Chem.MolFromSmarts("C(=O)O[#6]")
NITRILE_SMARTS = Chem.MolFromSmarts("C#N")
#allow F and Cl
HEAVY_HALOGEN_SMARTS = Chem.MolFromSmarts("[Br,I]")

#constants 
#number of examples that we want to keep per fragment (probably keep at 3)
MAX_EXAMPLES = 3

#not sure if we want to keep this or not, or find some way to ensure that we can still ultimately fully compose a molecule
#ensure fragment is not too big or too small
#especially on first pass, likely want to keep super large fragments, even if they end up basically becoming singletons
#fragments that are too small could likely just be treated as linkers
MIN_HEAVY_ATOMS = 5
MAX_HEAVY_ATOMS = 15

#definition of the empty fragment dictionary that we will populate
fragment_dict = {}

#likely remove, housekeeping from the prior workflow
#we will need to derive new handling for the 64b library
#this ocmmented block is how I iterated over the 2.6b library in a single superchunk
"""
#clobber and then create a folder by the superchunk name and move into it
os.system("rm -drf " + str(superchunk))
os.system("mkdir " + str(superchunk))

os.chdir(str(superchunk))

#iterate over all chunks and subchunks within the superchunk in the library to get the ligand fragments
#iterate chunks (0-99)
for i in range(0,100):

	#chunk base so add a preceeding 0 to the chunk if it is 0-9
	chunk_base = str(i)

	while len(chunk_base) < 2:
		chunk_base = "0" + chunk_base

	#derive chunk name, which is the superchunk followed by the chunk, with zeroes appended to the front to make a 5 digit string
	#i.e. superchunk 1 chunk 43 would be '00' '1' ''43' where 1 + 43 are concatenated, and then preceeding zeroes are appended
	chunk_id = str(superchunk) + chunk_base
	while len(chunk_id) < 5:
		chunk_id = "0" + chunk_id

	#adding a temporary test throttle so  Ican quickly check a a few chunks quickly
	#if chunk_id == "00008":
	#	break

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
"""

    #this is how I fed in a sdf file into rdkit, this will be different from what you do for feeding rdkit the smiles from the bz2
		supplier = Chem.SDMolSupplier("split_new_named_" + str(j) + ".sdf", removeHs=True)

    #enumerates over all ligands in the sdf file (will do some different type of iteration over bz2 smiles)
		for mol_idx, mol in enumerate(supplier):
        #safety net for if the current molecule is "none"
		    if mol is None:
		        continue

		    #derive the ligand name as an example
        #action will be different for bz2 since we have the name as a line entry, don't have to pull property from string b/c it doesn't exist
		    ligand_name = (
		        mol.GetProp("_Name")
		        if mol.HasProp("_Name")
		        else f"mol_{mol_idx}"
		    )

		    #adjust the ligand name to also include the superchunk, chunk, and subchunk that it came from
        #maybe could do something like this for the tracking for the up to 3 examples of the fragment/linker
		    ligand_name = str(superchunk) + ":" + chunk_id + ":" + str(j) + ":" + ligand_name


		    # Keep largest component (salt stripping)
        #realistically, this may not actually do anything, but is a good safety call
		    frags = Chem.GetMolFrags(mol, asMols=True)
		    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

        #molecular sanitization, keep; also probably doesn't actually do anything, but is a good safety net
		    Chem.SanitizeMol(mol)
		    Chem.Kekulize(mol, clearAromaticFlags=False)

		    # Find rotatable bonds
        #likely one of the most important things, we get the rotatable bonds based on the smarts the we define at the top
		    rot_bonds = mol.GetSubstructMatches(ROT_BOND_SMARTS)

        #makes an editable copy of the original molecule
		    em = Chem.EditableMol(mol)
        #iterates over the rotatable bonds in the molecule and removes them to get fragments
		    for a1, a2 in rot_bonds:
		        em.RemoveBond(a1, a2)

        #make a copy of the cut molecule
		    cut_mol = em.GetMol()
        #define rigid_frags as the fragments from the cut molecule
		    rigid_frags = Chem.GetMolFrags(
		        cut_mol, asMols=True, sanitizeFrags=True
		    )

        #set to indicate fragments that we hit and will skip if we have already seen the fragment
        #probably want to remove since we need to know every fragment that the ligand is composed of to compose a recipe; do not skip a found fragment, even if it is a duplicate
		    seen_this_ligand = set()

        #iterate over each fragment
		    for frag in rigid_frags:
		        # Size filter
            #we want to change this to have tiny fragments be established as a linker and just allow all huge fragments
		        if frag.GetNumHeavyAtoms() < MIN_HEAVY_ATOMS:
		            continue
            #comment
		        #if frag.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS:
		        #    continue

		        # Chemical exclusion filters
            #we do not want to blacklist problematic groups currently, so that we may compose the whole library
		        #if frag.HasSubstructMatch(ESTER_SMARTS):
		        #    continue
		        #if frag.HasSubstructMatch(NITRILE_SMARTS):
		        #    continue
		        #if frag.HasSubstructMatch(HEAVY_HALOGEN_SMARTS):
		        #    continue

		        #remove overly complex fragments with many rings or too many heteroatoms
            #also comment out, we do not want to remove overly complex fragments yet, just so we can account for the whole library
		        #num_rings = frag.GetRingInfo().NumRings()
		        #num_hetero = sum(1 for a in frag.GetAtoms() if a.GetAtomicNum() not in (1, 6))  # H + C only
		        #if num_rings > 4 or num_hetero > 8:
		        #    continue

            #derive the smiles of the fragment
		        frag_smiles = Chem.MolToSmiles(
		            frag,
		            canonical=True,
		            isomericSmiles=False
		        )

            #likely want to remove so that we retain all fragments
		        #if frag_smiles in seen_this_ligand:
		        #    continue
		        #seen_this_ligand.add(frag_smiles)

            #register a new fragment in the dictionart if it does not exist already and give the current ligand (with its source location) as an example
            #also initialize the count for the number of times the fragment was seen to 1
		        if frag_smiles not in fragment_dict:
		            fragment_dict[frag_smiles] = {
		                "count": 1,
		                "examples": [ligand_name]
		            }
		        #if the ligand was already in the dictionary, increment the encounter count and append the ligand example if we have not exceeded 3 (or whatever we use as the max example count)
            else:
		            fragment_dict[frag_smiles]["count"] += 1
		            if len(fragment_dict[frag_smiles]["examples"]) < MAX_EXAMPLES:
		                fragment_dict[frag_smiles]["examples"].append(ligand_name)
    #housekeeping cleanup from library processing
		#os.system("rm split_new_named_" + str(j) + ".sdf.tar.gz")
		#os.system("rm split_new_named_" + str(j) + ".sdf")



#loop over the dictionary once to remove odd fragments that are large and singleton fragments (1 occurenct and at least 10 heavy atoms)
#likely do not want to do this so that we can keep all fragments so that we may build all ligands from recipes
"""
for frag_smiles in list(fragment_dict.keys()):  # list() to allow deletion
	frag_info = fragment_dict[frag_smiles]

	if frag_info["count"] == 1:
		mol = Chem.MolFromSmiles(frag_smiles)
		if mol and mol.GetNumHeavyAtoms() >= 11:
			del fragment_dict[frag_smiles]
"""
#output the fragment dictionary to a csv file
#we will definitely want an output of the ligand fragment library
#will have to change things like file naming logic, but the overall concept holds
write_file = open(str(superchunk) + "_fragments.csv", "w")
write_file.write("fragment_smiles,occurences,example_ligands\n")
for frag in fragment_dict.keys():
	write_file.write(frag + "," + str(fragment_dict[frag]["count"]) + ",")
	#write the examples
	for ligname in fragment_dict[frag]["examples"]:
		write_file.write(str(ligname) + ";")

	write_file.write("\n")
