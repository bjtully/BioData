#!/usr/bin/python

"""Designed to take a full Genbank file (nucleotide scaffold with protein/CDS features) and convert it to a 
FASTA protein file. If there are multiple organisms in the Genbank file, multiple FASTA files will be generated.

The protein features must have either a locus_tag or accession number for the script to work.

Usage: python FullGenbank_to_MultipleFASTA.py <genbank.gbk>
"""

import sys
import string
from Bio import SeqIO

try:
	in_genbank = open(sys.argv[1], "r")
except IndexError:
	print "###Input full Genbank file to convert it to a FASTA protein file.\n\nUsage: python FullGenbank_to_MultipleFASTA.py <genbank.gbk>"


data_dict = {}

#process Genbank file
for record in SeqIO.parse(in_genbank, "genbank"):
	#Organism descript must be in the format ORGANISM NAME: ACCESSION NUMBER
	descp = record.description.split(":")[0]
	if "plasmid" in descp:
		continue
	seq_id = record.id
	#initial the dictionary. each k = organism accession number, v = a list with organism name at [0], and multiple dictionaries
	#of the putative CDS data at [1] to len(data_dict[seq_id])
	data_dict[seq_id] = [descp]
	for feature in record.features:
		#only interested in proteins
		if feature.type == "CDS":
			protein_dict = {}
			product_name = feature.qualifiers['product'][0]
			try:
				protein_seq = feature.qualifiers['translation'][0]
			except KeyError:
				continue
			#check to see if the genbank entry contains the protein accession ID
			if 'protein_id' in feature.qualifiers.keys():
				protein_id = feature.qualifiers['protein_id'][0]
				protein_dict[protein_id] = [product_name, protein_seq]
				#adds the entire putative CDS to the position in the list stored for each organism in data_dict
				data_dict[seq_id].append(protein_dict)
			#defaults to using the locus tag is accession is not present
			elif 'locus_tag' in feature.qualifiers.keys():
				locus_tag = feature.qualifiers['locus_tag'][0]
				protein_dict[locus_tag] = [product_name, protein_seq]
				#adds the entire putative CDS to the position in the list stored for each organism in data_dict
				data_dict[seq_id].append(protein_dict)
file_names = []			
#parse through each accession number in data_dict
for k in data_dict:
	#list contain the data for each genome
	genome = data_dict[k]
	#organism name stored in position [0]
	source_name = genome[0]
	#removes spaces in the organism name to create a filehandle name
	file_title = source_name.replace(" ", "")
	#Create a subset of all punctuation and use that to remove from organism names
	exclude = set(string.punctuation)
	file_title = ''.join(ch for ch in file_title if ch not in exclude)
	#creates a file for the organism using the name stored in description without spaces
	if file_title in file_names:
		output = open("%s.protein.faa" % file_title, "a")
	else:
		try:
			output = open("%s.protein.faa" % file_title, "w")
			file_names.append(file_title)	
		except IOError:
			continue
	#skips the organism name in position [0]. genome is a list of dictionaries. at each position starting at [1] until the
	#the length of the list
	for x in range(1,len(genome)):
		#dictionary at position [x]
		current_protein = genome[x]
		#for each protein create an entry in the FASTA file
		for s in current_protein:
			ident = s
			product_name = current_protein[s][0]
			protein_seq = current_protein[s][1]
			#ident can be either an accession ID or locus tag
			output.write(">%s %s|%s|%s\n%s\n" % (ident, product_name, source_name, k, protein_seq))
	output.close()
