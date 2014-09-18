Various Bioinformatic Tools
===========================

Tools that I have created to process bioinformatic data with potentially many purposes, not for a single project. If you are a member of the C-DEBI community, I can work with you to create generic scripts that process your data. 

##Full Genbank to Multiple FASTA

Script designed to process a GENBANK file that contains both the genomic scaffold sequence and the protein CDS features. 

Generally, this is the type of file output by IMG when multiple genomes are in the Genome Cart.

The purpose of this script is to generate a FASTA file of the proteins for each genome in the GBK file.

example.gbk contains a plasmid sequence in this format

CAUTION: The protein features MUST CONTAIN either a locus tag or accession number

###Dependencies

* [Biopython](http://biopython.org/wiki/Download)

###Usage
```
python FullGenbank_to_MultipleFASTA.py <genbank.gbk>
```
