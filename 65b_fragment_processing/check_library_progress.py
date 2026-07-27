#a simple script to parse over the library and determine which chunks are complete and which are not
#determined by checking the manifest and seeing if all batches are complete (have teh stats json file)
#will write a report on which chunks are complete and which are not

import os,sys

location = "/work/umassmed/thymelab_umassmed/65b_library_extracted_features/features/"

#open writers for complete and not complete
complete_writer = open(location + "complete_chunks.txt","w")
incomplete_writer = open(location + "incomplete_chunks.txt","w")

#write a header for the incomplete
incomplete_writer.write("chunk,total_batches_expected,batches_complete\n")

#iterate over the S and M folders
for r,d,f in os.walk(location):
	for dire in d:
		if dire.startswith("H") and ("M" in r or "S" in r):
			#we are at an individual directory, now look down the hcunk directory to determine its progress
			print(r + "/" + dire)

			#stats file counter
			stats_counter = 0

			#total expected batches from the manifest
			total_batches = 0

			for r2,d2,f2 in os.walk(r + "/" + dire):
				#get the manifest info and number of stats files to determine if there is a match
				for file in f2:


					if file == "manifest.tsv":
						manifest_reader = open(r2 + "/" + file, "r")
						
						for line in manifest_reader:
							#skip the header
							if line.startswith("task_id"):
								continue

							#get the first entry, which is the batch number; get the highest batch number listed
							line_batch = int(line.split()[0])
							if line_batch > total_batches:
								total_batches = line_batch

					#if the file is a batch file, increment the stats counter
					if file.startswith("batch_") and file.endswith(".stats.json"):
						stats_counter = stats_counter + 1

			#compare the expected batches with found complete batches
			if stats_counter == total_batches:
				complete_writer.write(r + "/" + dire + "\n")

			else:
				incomplete_writer.write(r + "/" + dire + "," + str(total_batches) + "," + str(stats_counter) + "\n")

