#!/usr/bin/python

"""The second iteration of a ARSC determining python script. This will be applied to test subset of 
genbank files from NCBI. It will determine the average ARSC for N, S, C, O and write it a text file.
This text file will be used going to foward to avoid calculating the same values multiple times
Further it will determine the ARSC for each protein in each genome, allowing for a future link
wherey a single point on the average ARSC - GC plot can be selected and each protein can be 
visualized."""

"""Must provide both the full genome and the sequences for each protein in FNA and GBK format"""

#THE TABLE for this analysis. Key = single letter amino acid code
#Value = numerical list with the order N - C - S - O - 1|2|3|4|5|6|7
#in position 5 a numerical code 1 = small, 2 = nucleophilic, 3 = hydrophobic, 4 = aromatic
#5 = acidic, 6 = amide, 7 = basic
amino_acid_atom_count = {"G":[0,0,0,0,1],
						 "A":[0,1,0,0,1], 
						 "S":[0,1,0,1,2],
						 "T":[0,2,0,1,2],
						 "C":[0,1,1,0,2],
						 "V":[0,3,0,0,3],
						 "L":[0,4,0,0,3],
						 "I":[0,4,0,0,3],
						 "M":[0,3,1,0,3],
						 "P":[0,3,0,0,3],
						 "F":[0,7,0,0,4],
						 "Y":[0,7,0,0,4],
						 "W":[1,9,0,0,4],
						 "D":[0,2,0,2,5],
						 "E":[0,3,0,2,5],
						 "N":[1,2,0,1,6],
						 "Q":[1,3,0,1,6],
						 "H":[2,4,0,0,7],
						 "K":[1,4,0,0,7],
						 "R":[3,4,0,0,7],						 	
}

from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab as pl
import sys
import os
import glob

#counts the elements in a protein a single amino acid at a time and returns a list of elements useage
def aa_count(protein):
	N = 0
	C = 0
	S = 0
	O = 0
	for x in range(0,len(protein)):
		current_amino_acid = protein[x]
		#checks for non-standard amino acid symbols
		if current_amino_acid in amino_acid_atom_count.keys():
			N += amino_acid_atom_count[current_amino_acid][0]
			C += amino_acid_atom_count[current_amino_acid][1]
			S += amino_acid_atom_count[current_amino_acid][2]
			O += amino_acid_atom_count[current_amino_acid][3]
	return [N, C, S, O]

#check to see if there is a file that contains ARSC calculations already
#if one does exist, pull out the Accession IDs to prevent the calculation from occurring twice
previous_calc_ids = []
#if not file, create one with the following tab delimited header line
if os.path.isfile("precalculated_avgARSC.txt") == False:
	quickout = open("precalculated_avgARSC.txt", "w")
	quickout.write("AccessionNumber	OrganismName	GC	AvgN_ARSC	AvgC_ARSC	AvgS_ARSC	AvgO_ARSC\n")
	quickout.close()
#if the file exist populate list to prevent repeat calculations
if os.path.isfile("precalculated_avgARSC.txt") == True:
	filehandle = open("precalculated_avgARSC.txt", "r")
	for line in filehandle:
		if line[:9] != "Accession":
			previous_calc_ids.append(line.split("\t")[0])
	filehandle.close()


#search the local directory and create 2 lists of files
#1 genbank file
#1 fasta format genomic scaffold file
gbk_names = glob.glob('*.gpff')
fna_names = glob.glob('*.genomic.fna')

#the major data storage dictionaries
protein_ARSC_dict = {}
genome_GC = {}

#parse through the list of genbank names
for n in gbk_names:
	#parse through each file
	for record in SeqIO.parse(open(n, "r"), "genbank"):
		#organisms Accession ID
		uniq_id = record.annotations['db_source'][18:-2]
		#check to see if it is in the file that has been previously precalculated
		if uniq_id not in previous_calc_ids:
			#sets organism name
			organism_name = record.annotations['source']
			#Protein Accessio number
			protein_id = record.id[:-2]
			protein_data = {}
			#counts the elements and places a list in a dictionary with Protein Accession as key
			protein_data[protein_id] = aa_count(record.seq)
			#appends the protein length to the list in the dictionary
			protein_data[protein_id].append(len(record.seq))
			#checks the validity of the dictionary
			if uniq_id not in protein_ARSC_dict.keys():
				#adds organism as a list if not present based on Organism Accession ID
				protein_ARSC_dict[uniq_id] = [organism_name]
				#appends the element values and protein length
				protein_ARSC_dict[uniq_id].append(protein_data)
			else:
				#if the Organism is already in the dictionary
				protein_ARSC_dict[uniq_id].append(protein_data)

#a list of all the organism Accession IDs from the ALL of the genbank files
active_list = protein_ARSC_dict.keys()
#parse through the list of FASTA files
for m in fna_names:
	#parse each FASTA file
	for record in SeqIO.parse(open(m, "r"), "fasta"):
		#pulls out the Accession ID
		a = record.id.split("|")
		fasta_id = a[3][:-2]
		#checks to see if the Accession was seen in the genbank files
		if fasta_id in active_list:
			#adds to dictionary with Organism Accession ID as key and GC content as value
			genome_GC[fasta_id] = GC(record.seq)

#append results to the output file
output = open("precalculated_avgARSC.txt", 'a')

#if not present create a file to store the individual protein data
if os.path.isfile("individual_protein_ARSC.txt") == False:
	quickout = open("individual_protein_ARSC.txt", "w")
	quickout.write("Protein Accession	Organism Accession	Organism Name	N ARSC	C ARSC	S ARSC	O ARSC	Protein Length\n")
	quickout.close()

#appends results to the protein output file
protein_output = open("individual_protein_ARSC.txt", "a")

#searching through the dictionary with all the data
for k in protein_ARSC_dict:
	#running counts
	total_residues = 0
	total_N = 0
	total_C = 0
	total_S = 0
	total_O = 0
	#k = Organism Accession ID
	#info_list = [Species Name, [N ARSC, C ARSC, S ARSC, O ARSC]]
	info_list = protein_ARSC_dict[k]
	#step through the list
	for x in info_list:
		#if the element is the species name, disregard it
		if isinstance(x, basestring) == True:
			continue
		#if it is the element information, step through and add values to the running counts
		else:
			for y in x:
				total_residues += x[y][4]
				total_N += x[y][0]
				total_C += x[y][1]
				total_S += x[y][2]
				total_O += x[y][3]
				#write out thes values for the protein_output file
				protein_output.write(str(y)+"\t"+str(k)+"\t"+str(protein_ARSC_dict[k][0])+"\t"+str(x[y][0])+"\t"+str(x[y][1])+"\t"+str(x[y][2])+"\t"+str(x[y][3])+"\t"+str(x[y][4])+"\n")
	#perform calculations for Avg ARSC values
	avg_N = float(total_N)/float(total_residues)
	avg_C = float(total_C)/float(total_residues)
	avg_S = float(total_S)/float(total_residues)
	avg_O = float(total_O)/float(total_residues)	
	#the corresponding %GC for the target genome
	gc = genome_GC[k]
	#the output of precalculated Avg. ARSC values
	output.write(str(k)+"\t"+str(protein_ARSC_dict[k][0])+"\t"+str(gc)+"\t"+str(avg_N)+"\t"+str(avg_C)+"\t"+str(avg_S)+"\t"+str(avg_O)+"\n")
















