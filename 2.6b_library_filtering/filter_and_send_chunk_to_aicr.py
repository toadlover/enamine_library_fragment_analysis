#the purpose of this script is to take a given chunk from the enamine library staged on the umass hpc, use the filtered csv report to remove duplicate enantiomers and very similar ligands, and then copy the chunk to the airc server

#imports
#this shouldn't need anything fancy like rdkit, just raw csv reads and file i/o and movement
#since this is a discreet and one-off operation, pathing will be hard-coded to make this fast to write

import os,sys

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
filtered_library_file = "/pi/summer.thyme-umw/2.6b_library_filtering/ari_test/" + working_superchunk + "/" + working_chunk + "/chunk_00000_similarity_report.csv.tar.gz"

#unfiltered library location
unfiltered_library_location = "/pi/summer.thyme-umw/enamine-REAL-2.6billion/" + working_superchunk + "/" + working_chunk

#copy the contents to the working location
os.system("cp -drf " + filtered_library_file + " " + working_location)
os.system("cp -drf " + unfiltered_library_location + " /condensed_params_and_db_*.tar.gz" + " " + working_location)

#unzip everything
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith(".gz"):
			os.system("tar -xzf " + file)

#go through the csv file and tidy it up

#remove shapedb data from each condensed folder and remove ligand params data that was flagged for removal for either being an enantiomer or too similar to another passed ligand

#compress everything

#send the files to aicr

#delete the files here for cleanliness