#!/usr/bin/python

"""Determining the file list within the ncbi refseq ftp using the following sets of scripts:
from ftplib import FTP
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login('anonymous', <EMAIL>)
ftp.cwd("/refseq/release/bacteria")
files = ftp.dir()
print files

Move this file structure to a new text file.

Parse text file for the sets of files -- i.e. only grab a files if both the *.genomic.fna.gz and
*.protein.faa.gz files are present.

Pull down target files."""

file_list = []
for line in open("/Users/btully/Documents/CDEBIPostdoc/DS-LA-03/directory_filelist.txt", "r"):
	a = line.split()
	file_list.append(a[8])

genomic_fna_list = []
protein_faa_list = []

for item in file_list:
	parts = item.split(".")
	if parts[1] == "nonredundant_protein":
		continue
	if "protein" and "gpff" in parts:
		protein_faa_list.append(item)
	if "genomic" and "fna" in parts:
		genomic_fna_list.append(item)

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

from ftplib import FTP
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login('anonymous', 'tully.bj@gmail.com')
ftp.cwd("/refseq/release/bacteria")

for file_id in download_list:
	testfile = open('%s' % file_id, 'wb' )
	ftp.retrbinary("RETR %s" % file_id, testfile.write )
	testfile.close()

