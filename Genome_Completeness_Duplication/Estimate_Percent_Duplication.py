#!/usr/bin/python
"""
Using the ouput of ConservedSingleCopyGenes.py -> singlecopygenes.faa this script will BLASTP
the proteins of an unknown organism to determine the (%) duplication within the genome,
good for SAGs and metagenome composite genomes.

A match to the CSCGs is determined using an 80(%) match alignment and 50(%) identity cutoff.

Usage: python Estimate_Percent_Duplication.py <CSCGs | singlecopygenes.faa> <target genome>.faa

Note: Especially if the (%) duplication is low, manual inspection will likely reveal that some or
all of the duplicate genes can be explained as the result of assembly errors/duplications
With manual inspection it may be possible to remove such errors from the predicted genome
"""

from Bio import SeqIO
import os
import sys


try:
	#CSCG file
	cscgs = str(sys.argv[1])
	#Target genome protein FASTA
	unknown = str(sys.argv[2])
except IndexError:
	print "Usage: python Estimate_Percent_Duplication.py <CSCGs | singlecopygenes.faa> <target genome>.faa"
	exit()


#Runs makeblastdb if this is the first time singlecopygene.faa is used as
#BLAST database
#If a new version of singlecopygenes.faa is used be sure to remove old
#BLAST indices
if os.path.isfile("singlecopygenes.faa.phr") == False:
	os.system("makeblastdb -in %s -dbtype prot" % cscgs)
#Run a BLASTP with e = 0.001, 1 best match, & tab output with the designated
#values
os.system("blastp -query %s -db %s -out BLAST_%s_%s -evalue 0.001 -max_target_seqs 1 \
		-outfmt '6 qseqid qlen qstart qend sseqid stitle slen sstart send bitscore pident evalue'" % (unknown, cscgs, cscgs, unknown))
 
#Simple count of the number of CSCGs in the protein FASTA file
num_cscgs = 0
for record in SeqIO.parse(open(cscgs, "r"), "fasta"):
	num_cscgs +=1

#a version of the function from ConservedSingleCopyGenes.py
def determine_ifCSCG(data):
	#Percent length of the subject in the CSCG BLAST database
	perclength = int((float(data[8])-float(data[7])+1)/float(data[6])*100)
	pident = float(data[10])
	if int(perclength) >= 80 and int(pident) >=50:
		return data[4]

#open the BLAST output
file_handle = open("BLAST_%s_%s" % (cscgs, unknown), "r")
all_identified_cscgs = []
#for each line in the BLAST output, determine if the protein is a match to a CSCG
#this results in a list of the singlecopygenes.faa sequence names that had matches
#in the target organism
for line in file_handle:
	a = line.split("\t")
	if determine_ifCSCG(a) != None:
		all_identified_cscgs.append(determine_ifCSCG(a))
		#print line
		#un-comment the above line to see the matches - can be useful for
		#determining if the duplication is real or assembly artifact

dupe_count = 0
#Use set to remove all the duplicates in the list
#Parse through the set of protein ids
for i in set(all_identified_cscgs):
	#then if the protein id occurs more than once in the list with
	#duplicated allowed - increase the duplicate count
	if all_identified_cscgs.count(i) > 1:
		print i
		dupe_count += 1

#Final duplication estimate determined by number of identified duplicates
#divided by the total number of CSCGs
perc_dupe = float(dupe_count)/float(num_cscgs)*100
print "Percent duplication is: " + str(perc_dupe)