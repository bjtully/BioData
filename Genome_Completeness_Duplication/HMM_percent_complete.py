#!/usr/bin/python

"""
coregenome_tigrfam.txt generated using HMM_tbl_parse.py provides the expected number of conserved TIGRFAM functions for
a particular group of related organisms.

Using HMMER3, a genome of unknown completeness is compared to coregenome_tigrfam.txt to determine the (%) completeness 

HMMER paramters:
hmmscan -o <name>.out --cpu 25 --noali -E 0.00001 --tblout <name>.tbl /location/of/TIGRFAMs_14.0_HMM.LIB 
	<input unknown PROTEIN FAA>

Usage: python HMM_percent_complete.py coregenome_tigrfam.txt <target>.tbl
"""

import sys
from collections import Counter

try:
	coregenome = open(str(sys.argv[1]), "r")
	unknown_tbl = open(str(sys.argv[2]), "r")
except IndexError:
	print "Usage: python HMM_percent_complete.py coregenome_tigrfam.txt <target>.tbl"
	exit()

core = {}
expected_total = 0

#counts the number of TIGRFAM occurences in coregenome_tigrfam.txt
#populates dictionary that will be used to compare to the results of the
#unknown organism
for line in coregenome:
	a = line.split("\t")
	core[a[0]] = a[1]
	expected_total += int(a[1])

genome_storage = {}


#similar to HMM_tbl_parse.py
#Find the best TIGRFAM match for each protein in the unknown genome
for line in unknown_tbl:
	if line[0] != "#":
		a = line.split()
		tigr = a[0]
		query_name = a[2]
		bit_score = float(a[5])
#Compare bit scores to determine if the current match is the best match
		if query_name in genome_storage.keys():
			if bit_score > genome_storage[query_name][0]:
				genome_storage[query_name] = [bit_score]
				genome_storage[query_name].append(tigr)
			else:
				continue
		if query_name not in genome_storage.keys():
			genome_storage[query_name] = [bit_score]
			genome_storage[query_name].append(tigr)
#genome_storage[protein_id] = [bit_score, TIGRFAM]

#genome_tigr = a list of all the match TIGRFAMs 
genome_tigr = []
for k in genome_storage:
	genome_tigr.append(genome_storage[k][1])

#A dictionary representing the counts of each TIGRFAM for the unknown genome
genome_cnt = Counter(genome_tigr)

#When the unknown genome matches a TIGRFAM of the core genome, determine
#the minimum number for that TIGRFAM
#Ex. if unknown = 5, but core = 2, add 2 to the count b/c core genome = 2
#Ex. if unknown = 2, but core = 5, add 2 to the count b/c unknown is missing 3
identified_total = 0
for k in genome_cnt:
	if k in core.keys():
		if int(genome_cnt[k]) > int(core[k]):
			identified_total += int(core[k])
		if int(genome_cnt[k]) <= int(core[k]):
			identified_total += int(genome_cnt[k])

#Simple calculation
perc_complete = float(identified_total)/float(expected_total)*100

print "Percent complete: " + str(perc_complete)
