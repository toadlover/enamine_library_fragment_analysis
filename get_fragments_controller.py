#the purpose of this script is to fire off the script to get fragments from a superchunk for all 530 superchunks, and operate from the LSF job distributor
#imports
import os,sys

for i in range(0,531):
	print(i)
	os.system("sleep 1")
	command = "bsub -q short -W 8:00 -u \"\" -R \"rusage[mem=10000]\" \"python get_fragments_from_superchunks.py " + str(i) + "\""
	print(command)
	#os.system(command)