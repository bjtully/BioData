#!/usr/bin/python

"""
Usage: python ftpgrab_refseq.py <YOUR EMAIL> <bacteria or archaea>
Pull down target protein GENBANK and genomic FASTA files."""

import sys

#Access NCBI FTP. Login. Change Refseq directory to either bacteria of archaea. Get of list of the files in the directory.
from ftplib import FTP
ftp = FTP('ftp.ncbi.nih.gov')
try:
	ftp.login('anonymous', "%s" % str(sys.argv[1]))
except IndexError:
	print "Please provide correct input.\nUsage: python ftpgrab_refseq.py <YOUR EMAIL> <bacteria or archaea>"
	exit()
try:
	ftp.cwd("/refseq/release/%s" % str(sys.argv[2]).lower())
except IndexError:
	print "Please provide correct input.\nUsage: python ftpgrab_refseq.py <YOUR EMAIL> <bacteria or archaea>"
	exit()
raw_file_list = []
ftp.dir(raw_file_list.append)

genomic_fna_list = []
protein_faa_list = []

#Convert directory information into file names
file_list = []
for line in raw_file_list:
	a = line.split()
	file_list.append(a[8])

#Check each file name for protein genbank and genomic fasta status.
for item in file_list:
	parts = item.split(".")
	if parts[1] == "nonredundant_protein":
		continue
	if "protein" and "gpff" in parts:
		protein_faa_list.append(item)
	if "genomic" and "fna" in parts:
		genomic_fna_list.append(item)

#Create a list of all of the genomic. Check if the protein file is in the genomic list. Using the protein list
#come up with a download list for the genomic files.
genomic_ids = []
protein_ids = []
download_list = []
for name in genomic_fna_list:
	b = name.split(".")
	genomic_ids.append(b[1])
for name in protein_faa_list:
	b = name.split(".")
	if b[1] in genomic_ids:
		protein_ids.append(b[1])
		download_list.append(name)
for name in genomic_fna_list:
	b = name.split(".")
	if b[1] in protein_ids:
		download_list.append(name)

#Connect to FTP
from ftplib import FTP
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login('anonymous', '%s' % str(sys.argv[1]))
ftp.cwd("/refseq/release/%s" % str(sys.argv[2]).lower())

#Download target files.
for file_id in download_list:
	testfile = open('%s' % file_id, 'wb' )
	ftp.retrbinary("RETR %s" % file_id, testfile.write )
	testfile.close()