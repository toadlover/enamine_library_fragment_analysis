#a basic script written with help from chatgpt to test a basic operation for pulling fragments from a single test sdf file in the 2.6B ligand library
from rdkit import Chem
from collections import defaultdict

ROT_BOND_SMARTS = Chem.MolFromSmarts(
    "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"
)

ESTER_SMARTS = Chem.MolFromSmarts("C(=O)O[#6]")
NITRILE_SMARTS = Chem.MolFromSmarts("C#N")
HEAVY_HALOGEN_SMARTS = Chem.MolFromSmarts("[Cl,Br,I]")

MAX_EXAMPLES = 3
MIN_HEAVY_ATOMS = 3

fragment_dict = {}

supplier = Chem.SDMolSupplier("test.sdf", removeHs=True)

for mol_idx, mol in enumerate(supplier):
    if mol is None:
        continue

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
                "examples": [mol_idx]
            }
        else:
            fragment_dict[frag_smiles]["count"] += 1
            if len(fragment_dict[frag_smiles]["examples"]) < MAX_EXAMPLES:
                fragment_dict[frag_smiles]["examples"].append(mol_idx)


#output the fragment dictionary to a csv file
write_file = open("test_file_fragments.csv", "w")
write_file.write("fragment_smiles,occurences,example_ligands\n")
for frag in fragment_dict.keys():
	write_file.write(frag + "," + str(fragment_dict[frag]["count"]) + ",")
	#write the examples
	for ligname in fragment_dict[frag]["examples"]:
		write_file.write(str(ligname) + ";")

	write_file.write("\n")
