#!/usr/bin/python

"""
Randomly subsample an interleaved FASTQ file to generate a smaller sample set.

In order to save time, this script can be used to make mutliple subsets of
different read numbers from one execution. However, only a single shuffle 
event occurs for each pass through the script. For example, if you want 100
and 200 reads sampled from 1,000 possible sequences, this script will randomize
all of the possible sequences. The 100 subsample will be random. The 200 
subsample will be identical over the first 100 reads, but the second 100 reads
will be random.

If you want two subsamples of the same number, but each to be a different subset
of random reads, the script needs to be run twice.

All reads are read into memory, so for large datasets be sure to have sufficient
memory.

Usage: python FASTA_subsample.py <input FASTA> <list of increasing desired subsample cutoffs>
Ex: python FASTA_subsample.py test.fasta 1000 2000 4000 8000
Will produce 4 interleaved, subsample FASTA files with the defined number of
reads.
"""

from Bio import SeqIO
from random import shuffle
import sys,os

#Set input name and open target file
try:
	file_name = str(sys.argv[1])
	filehandle = open(file_name, "r")
except IndexError:
	print "Usage: python FASTA_subsample.py <input FASTA> <list of increasing desired subsample cutoffs>"
	exit()
#A list of at least 1 subsample length
number_subsample = sys.argv[2:]


#The part works based on the assumption that the FASTA is correctly interleaved
#We avoid checking the key list of total_seq_list by saving the sequence ID
#and using the intervleaved nature of the file to produce the dictionary faster
total_seq_list = {}
previous_id = ""
for record in SeqIO.parse(filehandle, "fasta"):
	#If you are looking at the second sequence of the pair, append parts
	if record.id == previous_id:
		total_seq_list[record.id].append(record.description)
		total_seq_list[record.id].append(record.seq)
	#If this is a new sequence, start a new list in the dictionary
	#set the holder value to the correct value
	if record.id != previous_id:
		total_seq_list[record.id] = []
		total_seq_list[record.id].append(record.description)
		total_seq_list[record.id].append(record.seq)
		previous_id = record.id
filehandle.close()

#Perform the action of popluation at keys list once
all_ids = total_seq_list.keys()
#Output feedback
print "Number of sequence pairs: " + str(len(all_ids))

#Perform the randomization of the sequences
#Shuffles the list of all the sequence IDs to provide a new order of sequences
shuffle(all_ids)

#For each target number of reads defined by the user
for x in range(len(number_subsample)):
	#This portion of the script is based on some assumptions
	#1. Final output FASTA files will have the number of the subsample appended
	#to the end of the file name
	#2. The first subsample will be created
	#3. After the first subsample, each additional subsample will use the previous
	#file to create the next larger subsample
	
	#Check to see if the first subsample FASTA exists
	if os.path.isfile("%s_%s" % (file_name, number_subsample[0])) == True:
		#Create the new output file
		output = open("%s_%s" % (file_name, number_subsample[x]), "w")
		#Read in and write all of the lines from the previous subsample FASTA
		for line in open("%s_%s" %(file_name, number_subsample[x-1]), "r"):
			output.write(str(line))
		#Determine the which sequences need to be added to the output file
		#based on where the previous file ends and the number of total needed
		pairs_start = int(number_subsample[x-1])/2
		pairs_end = int(number_subsample[x])/2
		#Using the shuffled list of sequence IDs, read through until target
		#number is reached, accounting for overlaps in numbers
		for target_seq in all_ids[pairs_start+x:pairs_end+x]:
			#Read data for each sequence pair
			target_data = total_seq_list[target_seq]
			seq_name1 = target_data[0]
			seq1 = target_data[1]
			seq_name2 = target_data[2]
			seq2 = target_data[3]
			#Printed to outfile
			output.write(">"+str(seq_name1)+"\n"+str(seq1)+"\n>"+str(seq_name2)+"\n"+
				str(seq2)+"\n")
		output.close()
		print "File %s_%s created." % (file_name, number_subsample[x])

	else:	
		#Creates first output file
		output = open("%s_%s" % (file_name, number_subsample[x]), "w")
		#Determine how many pairs of reads are required
		#Target number of sequences divided by 2
		num_pairs = int(number_subsample[x])/2
		#Using the shuffled list of sequence IDs, read through until target
		#number is reached
		for target_seq in all_ids[:num_pairs]:
			#Read data for each sequence pair
			target_data = total_seq_list[target_seq]
			seq_name1 = target_data[0]
			seq1 = target_data[1]
			seq_name2 = target_data[2]
			seq2 = target_data[3]
			#Printed to outfile
			output.write(">"+str(seq_name1)+"\n"+str(seq1)+"\n>"+str(seq_name2)+"\n"+
				str(seq2)+"\n")
		output.close()
		print "File %s_%s created." % (file_name, number_subsample[x])
