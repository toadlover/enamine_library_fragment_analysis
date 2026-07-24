#the purpose of this script is to take a given chunk from the enamine library staged on the umass hpc, use the filtered csv report to remove duplicate enantiomers and very similar ligands, and then copy the chunk to the airc server

#imports
#since this is a discreet and one-off operation, pathing will be hard-coded to make this fast to write

import os,sys

from rdkit import Chem
from rdkit.Chem import AllChem

#only arg needed is the chunk, and the processing can go from there
working_chunk = str(sys.argv[1])

#if the working chunk does not have enough leading zeroes (should be a 5 digit value), append leading zeroes
while len(working_chunk) < 5:
	working_chunk = "0" + working_chunk

#derive the working superchunk; cast to int to quickly cut leading zeroes and then cast back to str
working_superchunk = str(int(str(working_chunk[0:3])))

print("superchunk:", working_superchunk)
print("chunk:", working_chunk)

#working location
working_location = "/pi/summer.thyme-umw/2.6b_library_filtering/filtering_space/" + working_superchunk + "/" + working_chunk

#make the location
os.system("mkdir -p " + working_location)

#os.system("sleep 1")

#move here
os.chdir(working_location)

print(os.getcwd())

#filtered library file location
filtered_library_file = "/pi/summer.thyme-umw/2.6b_library_filtering/ari_test/" + working_superchunk + "/" + working_chunk + "/chunk_" + working_chunk + "_similarity_report.csv.tar.gz"

#unfiltered library location
unfiltered_library_location = "/pi/summer.thyme-umw/enamine-REAL-2.6billion/" + working_superchunk + "/" + working_chunk

#copy the contents to the working location
print("copying files")
os.system("cp -drf " + filtered_library_file + " " + working_location)
os.system("cp -drf " + unfiltered_library_location + "/condensed_params_and_db_*.tar.gz" + " " + working_location)
#os.system("cp -drf " + unfiltered_library_location + "/split_new_named_*.tar.gz" + " " + working_location)
os.system("cp -drf " + unfiltered_library_location + "/blacklist_file.csv" + " " + working_location)
print("unzipping files")
#unzip everything
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith(".gz"):
			os.system("tar -xzf " + file)

#delete the copied tar files since we do not need the originals
os.system("rm *gz")

#go through the csv files and tidy it up, write to a new file; merge the blacklist filtering with the duplicate filtering
print("identifying blacklist ligands")
blacklist_ligands = []
read_file = open("blacklist_file.csv","r")
for line in read_file.readlines():
	lig = line.split(",")[0]
	blacklist_ligands.append(lig)
	#if a ligand is in the blacklist, delete its data from the condensed params
	os.system("rm condensed*/single_conf_params/" + lig + "_shorthand_params.txt")
read_file.close()

#now, read through the report file
print("removing ligands, writing report, and sdf files")

read_file = open("chunk_" + working_chunk + "_similarity_report.csv","r")
#actually, we probably don't even need a manifest to write. just make the sdfs
#write_file = open("chunk_" + working_chunk + "_manifest.csv","w")

#make folders for where to host the shorthand params and sdfs
os.system("mkdir shorthand_params sdfs")

#read the report file and process for ligand keeping 
for line in read_file.readlines():
	#skip the header
	if line.startswith("global_"):
		continue

	#get the ligand name and smiles string
	lig = line.split(",")[1]
	smiles = line.split(",")[11]
	retained = line.split(",")[2]

	#if the ligand is in the blacklist, continue
	if lig in blacklist_ligands:
		continue

	#if the ligand was not retained, continue
	if retained == "removed":
		continue

	#otherwise, move the condensed params and make an sdf
	os.system("mv condensed*/single_conf_params/" + lig + "_shorthand_params.txt shorthand_params")

	mol = Chem.MolFromSmiles(smiles)
	if mol is None:
		raise ValueError(f"Invalid SMILES string: {smiles}")

	# Add explicit hydrogens before generating a 3D conformer.
	mol = Chem.AddHs(mol)

	# Store the ligand name in the SDF molecule-title field.
	mol.SetProp("_Name", lig)

	# Generate and optimize a 3D conformer.
	embed_status = AllChem.EmbedMolecule(
		mol,
		randomSeed=42,
		useRandomCoords=True,
	)
	if embed_status != 0:
		raise RuntimeError(f"Could not generate a 3D conformer for {lig}")

	# MMFF is preferred when parameters are available; otherwise use UFF.
	if AllChem.MMFFHasAllMoleculeParams(mol):
		AllChem.MMFFOptimizeMolecule(mol)
	else:
		AllChem.UFFOptimizeMolecule(mol)

	# Write the molecule to an SDF file.
	writer = Chem.SDWriter("sdfs/" + lig + ".sdf")
	writer.write(mol)
	writer.close()	



#remove the condensed folder and csvs
print("condensed cleaning")
os.system("rm -drf condensed* *.csv*")

#compress everything
print("compressing data to transfer")
os.system("tar -czf sdfs.tar.gz sdfs")
os.system("tar -czf shorthand_params.tar.gz shorthand_params")

#send the files to aicr
print("sending to aicr")
#"/work/umassmed/thymelab_umassmed/2.6b_conformer_library/" + working_superchunk + "/" + working_chunk
os.system("ssh -i ~/.ssh/id_ed25519_aicr ari_ginsparg_umassmed@login.aicr.ai 'mkdir -p /work/umassmed/thymelab_umassmed/2.6b_conformer_library/" + working_superchunk + "/" + working_chunk + "/' && scp -i ~/.ssh/id_ed25519_aicr sdfs.tar.gz shorthand_params.tar.gz ari_ginsparg_umassmed@login.aicr.ai:/work/umassmed/thymelab_umassmed/2.6b_conformer_library/" + working_superchunk + "/" + working_chunk)
#delete the files here for cleanliness
#os.system("rm -drf *")