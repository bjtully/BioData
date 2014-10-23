Various Bioinformatic Tools
===========================

Tools that I have created to process bioinformatic data with uses beyond a single project. 
If you are a member of the C-DEBI community, I can work with you to create generic scripts that process your data. 

##Full Genbank to Multiple FASTA

Script designed to process a GENBANK file that contains both the genomic scaffold sequence and the protein CDS features. 

Generally, this is the type of file output by IMG when multiple genomes are in the Genome Cart.

The purpose of this script is to generate a FASTA file of the proteins for each genome in the GBK file.

example.gbk contains an example sequence in this Genbank format.

CAUTION: The protein features MUST CONTAIN either a locus tag or accession number

###Dependencies

* [Biopython](http://biopython.org/wiki/Download)

###Usage
```
python FullGenbank_to_MultipleFASTA.py <genbank.gbk>
```

##Random Subsample of FASTA File

Script designed to generate random subsamples of a given FASTA file.

Assumptions of the script:

* Sequences are paired-end reads

* FASTA file is interleaved

The script allows to the generation of multiple subsamples from the same FASTA file. For example, if you want FASTA files with subsamples of 100, 200, 400, and 800 reads from a FASTA file with 1,500 reads, this script can produce this in a single pass.

To increase the speed of the script, subsequent subsamples are not completely random relative to the other subsamples. In the example above. Subsample-100 is a FASTA that contains a random subset of the 1,500 reads. In the FASTA for subsample-200, the first 100 sequences match the sequences in the subsample-100 FASTA, but the second 100 sequences are a random set. So on, for subsample-400 and subsample-800.

If you want to random set of reads from the same FASTA, run the script multiple times.

CAUTION: All of the sequences are read into memory, so for large datasets plan accordingly.

###Dependencies

* [Biopython](http://biopython.org/wiki/Download)

###Usage
```
python FASTA_subsample.py <input FASTA> <list of increasing desired subsample cutoffs>
EX. python FASTA_subsample.py test.fasta 100 200 400 800 
```