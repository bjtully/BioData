#!/usr/bin/python

"""
Determine the conserved single copy genes for a group of genomes. It is recommended to determine which genome represent the
"most central" genome of the group of interest. This can be be done by creating a phylogenetic tree using a marker gene
from each genome.

An iterative BLASTP starting with the user determined first genome will be used to to determine single copy genes (occuring 
only once per genome) with at least 80(%) alignment length and a user determined (%) AAID. A final protein FASTA file is 
is generated that then can be used with Estimate_Percent_Duplication.py to determine (%) duplication.

Usage: python ConservedSingleCopyGenes.py <phylogenetically central genome>.faa <percent AAID cutoff>


"""

from Bio import SeqIO
import glob
import os
from random import shuffle
import sys

#Read in all of the protein.faa files in the directory
#These file can be generated from a Genbank with FullGenbank_to_MultipleFASTA.py
proteome_list = glob.glob("*.protein.faa")
if len(proteome_list) == 0:
	print "Please execute ConservedSingleCopyGenes.py in a directory with the target protein.faa files"
	exit()
#The user defined genome is self-BLAST against itself to determine the initial list of single copy genes
try:
	first_genome = str(sys.argv[1])
	cutoff = sys.argv[2]
except IndexError:
	print "Usage: python ConservedSingleCopyGenes.py <phylogenetically central genome>.faa <percent AAID cutoff>"
	exit()

#Remove that first genome from the list of genomes that will be iteratively processed
proteome_list.remove(first_genome)

#function that can determine duplicates or matches, depending on the target
#data = a list from the BLAST output below
#check = if the paritcular call is looking for duplicates or matches
#value = user defined (%)AAID
def determine_ifDuplicate(data, check, value):
	if str(data[0]) != str(data[4]):
		if str(check) == "dupes":
			#percent alignment length
			perclength = int((float(data[8])-float(data[7])+1)/float(data[6])*100)
			#percent AAID
			pident = float(data[10])
			#the self BLAST is limited to matches with 80 (%) alignment length
			#and 80 (%) AAID
			if int(perclength) >= 80.0 and int(pident) >= 80.0:
				return data[4]
		if str(check) == "match":
			perclength = int((float(data[3])-float(data[2])+1)/float(data[1])*100)
			pident = float(data[10])
			#the iterative BLAST is limited to matches with 80 (%) alignment length
			#and user defined (%) AAID
			if int(perclength) >= 80.0 and int(pident) >= int(value):
				return data[0]

#reads in a FASTA sequence. In this case, the singlecopygenes.faa file created after each
#successive BLAST
#file_name = in this case is singlecopygenes.faa
#duplications = a list of duplications determined by the BLAST output
def remove_Dupes(file_name, duplications):
	records = []
	for record in SeqIO.parse(open(file_name, "r"), "fasta"):
		if str(record.id) not in duplications:
			records.append(record)
	#records of a generally smaller singlecopygenes.faa file
	return records

#reads in a FASTA sequence. In this case, the singlecopygenes.faa file created after each
#successive BLAST
#file_name = in this case is singlecopygenes.faa
#matches = a list of matches determined by the BLAST output
def keep_Matches(file_name, matches):
	records = []
	for record in SeqIO.parse(open(file_name, "r"), "fasta"):
		if str(record.id) in matches:
			records.append(record)
	return records

#if singlecopygenes.faa has never been generated, this will perform the self BLAST
#against the user defined first genome
if os.path.isfile("singlecopygenes.faa") == False:
	#makes BLAST directory
	os.system("makeblastdb -in %s -dbtype prot" % first_genome)
	#performs BLASTP of the first genome against itself
	#the out format is is tab format with the specific values listed in the results
	os.system("blastp -query %s -db %s -out BLAST_methano1_self -evalue 0.001 -max_target_seqs 5 \
		-outfmt '6 qseqid qlen qstart qend sseqid stitle slen sstart send bitscore pident evalue'" % (first_genome, first_genome))

	dup_list = []
	#opens the self BLAST output
	file_handle = open("BLAST_methano1_self", "r")
	#splits the line based on tab and process each list of the line with the function determine_ifDuplicate
	for line in file_handle:
		#input list
		a = line.split("\t")
		#the third input to function is fixed as the value is set at 80(%)
		#checks to see if the result is positive
		if determine_ifDuplicate(a, 'dupes', 0) != None:
			#adds returned values to duplication list
			dup_list.append(determine_ifDuplicate(a, 'dupes', 0))
	file_handle.close()
	#creates record file to create singlecopygenes.faa using fucntion remove_Dupes
	out_records = remove_Dupes(first_genome, dup_list)
	#create FASTA file
	SeqIO.write(out_records, "singlecopygenes.faa", "fasta")

#once singlecopygenes.faa is created
if os.path.isfile("singlecopygenes.faa") == True:
	#process through each genome in the genome list
	for x in range(0,len(proteome_list)):
		current_genome = proteome_list[x]
		#current number of the genome is the current index+1
		current_number = x+1
		#make BLAST database
		os.system("makeblastdb -in %s -dbtype prot" % current_genome)
		#Perform same BLAST as above
		os.system("blastp -query singlecopygenes.faa -db %s -out BLAST_singlecopy_methano%s -evalue 0.001 -max_target_seqs 5 \
			-outfmt '6 qseqid qlen qstart qend sseqid stitle slen sstart send bitscore pident evalue'" % (current_genome, current_number))
		match_list = []
		dup_list = []
		#parse through BLAST output
		file_handle = open("BLAST_singlecopy_methano%s" % current_number, "r")
		for line in file_handle:
			a = line.split("\t")
			#looking for matches. passes along the user defined (%) AAID
			if determine_ifDuplicate(a, 'match', int(cutoff)) != None:
				#checks to see is the return match already occurs in the matched list
				#if so, this indicates that it is not single copy in the current genome
				#and needs to be excluded from the next version of singlecopygenes.py
				if determine_ifDuplicate(a, 'match', int(cutoff)) in match_list:
					dup_list.append(determine_ifDuplicate(a, 'match', int(cutoff)))
				if determine_ifDuplicate(a, 'match', int(cutoff)) not in match_list:
					match_list.append(determine_ifDuplicate(a, 'match', int(cutoff)))
		file_handle.close()
		#removes sequences in the dup_list from match_list
		for s in dup_list:
			if s in match_list:
				match_list.remove(s)
		#creates record file to the new version of singlecopygenes.faa using keep_Matches
		out_records = keep_Matches("singlecopygenes.faa", match_list)
		SeqIO.write(out_records, "singlecopygenes.faa", "fasta")
		#prints a statement of current progress
		print "Updated singlecopygenes.faa", current_genome, len(match_list)


