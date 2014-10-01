#!/usr/bin/python

"""Parse through the TBL output format of HMMER3 comparison between TIGRFAM and proteome of multiple genomes.
The output will be a list of TIGRFAMs that represent the 'core' genome of the group and their representative
counts.

The inpur for this script is all HMMER3 results in TBL format in the desired directory.

A list is made which contains the minimum number of TIGRFAM hits common for all of the input genomes.

For example: if genome 1 has 3 copies of the TIGRFAM00001 and genome 2 has 2 copes of TIGRFAM00001, then the common shared 
value would be 2 for the core genome.

Produces an fixed output called: coregenome_tigrfam.txt

Collection of HMMER3 outputs, produced using the following parameters:
hmmscan -o <name>.out --cpu 25 --noali -E 0.00001 --tblout <name>.tbl /location/of/TIGRFAMs_14.0_HMM.LIB <input PROTEIN FAA>

Usage: python HMM_tbl_parse.py
test
"""

import glob
from collections import Counter

#Collect all tables
tbl_names = glob.glob('*.tbl')
#simple check to see if the current directory contains any TBL outputs
if len(tbl_names) == 0:
	print """Please execute HMM_tbl_parse.py in a directory with TBL output(s) from HMMER3.
TBL outputs can be produced using HMMER3:
hmmscan -o <name>.out --cpu 25 --noali -E 0.00001 --tblout <name>.tbl /location/of/TIGRFAMs_14.0_HMM.LIB <input PROTEIN FAA>"""
	exit()

data_storage = {}

#Parse through the names of tables
for tbl_in in tbl_names:
	filehandle = open(tbl_in, "r")
	#The best TIGRFAM match for each proteon for each genome stored in this dictionary
	genome_storage = {}
	for line in filehandle:
		if line[0] != "#":
			a = line.split()
			#Hit TIGRFAM
			tigr = a[0]
			#Unique query name of the input genome protein
			#Using the input FAA sequences from FullGenbank_to_MultipleFASTA.py
			#this ID is the Accession number
			query_name = a[2]
			#Used to compare TIGRFAM matches to ensure the best match is stored
			bit_score = float(a[5])
			#if the protein is a already in the dictionary, ensure that the best  
			#hit determined by bit score is stored for proteins that hit 
			#multiple TIGRFAMs
			if query_name in genome_storage.keys():
				if bit_score > genome_storage[query_name][0]:
					genome_storage[query_name] = [bit_score]
					genome_storage[query_name].append(tigr)
#genome_storage[protein_id] = [bit_score, TIGRFAM]
				else:
					continue
			#if protein not in dictionary, key = protein query_id
			#value = list of bitscore and matching TIGRFAM
			if query_name not in genome_storage.keys():
				genome_storage[query_name] = [bit_score]
				genome_storage[query_name].append(tigr)
	#Parse the input TBL name at the <name>.faa
	#the unique species name in the file title is used as a uniq key in overall
	#dictionary
	org_id = tbl_in.split(".")[0]
	data_storage[org_id] = []
	#each organisms represented as a input TBL file is added to as a key to the
	#dictionary, data_storage, with values = list of matching TIGRFAMs
	for k in genome_storage:
		data_storage[org_id].append(genome_storage[k][1])
#data_storage[organismname] = [TIGRFAM00001, TIGRFAM00002, TIGRFAM00001, ...]


#stores the results of the minimum number of TIGRFAM occurrences
master_list = {}


for k in data_storage:
	#z = list of TIGRFAM matches for a particular genome
	z = data_storage[k]
	#Counter creates a dictionary where the keys = TIGRFAM IDs and
	#values = the number of occurrences for that ID
	data_list = Counter(z)
	new_master = {}
	#if master_list is empty, populate it with the first genome
	if bool(master_list) == False:
		master_list = data_list
	else:
	#TIGRFAMs have to be in all genomes, so only concerned with matches in for
	#each genome
		for x in data_list:
			#check to see if TIGRFAM ID is in master list
			if x in master_list.keys():
				#find the lowest number of occurences for the TIGRFAM
				#add to new_master, which will replace the master_list
				#after each genome is assessed 
				if int(data_list[x]) <= int(master_list[x]):
					new_master[x] = data_list[x]
				if int(data_list[x]) > int(master_list[x]):
					new_master[x] = master_list[x]
		master_list = new_master

output = open('coregenome_tigrfam.txt', "w")
#print the final ouput as TIGRFAM \t minimum number of occurrences
for k in master_list:
	output.write(str(k)+"\t"+str(master_list[k])+"\n")
output.close()





	









