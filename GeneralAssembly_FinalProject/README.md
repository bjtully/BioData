Getting Started
===============

An on going data science project

##Bacteria and/or Archaea Files
The script ftpgrab_refseq.py can grab all of the protein GENBANK and corresponding genomic scaffold FASTA files from the NCBI RefSeq database. CAUTION: Be sure to have 50-100 GB of space free.

###Usage
```
python ftpgrab_refseq.py <YOUR@EMAILADDRESS.com> <bacteria or archaea>
```

##Multi Genome ARSC
multiGenomeARSC_V2.py is the current iteration of a script that will determine the occurrence of elements (N, C, O, S) in the side chains of proteins for a bacterial/archaeal genome.

The counting method relies on counts of atoms per residue side chain (ARSC) (1,2). An average value is determined for each element by counting each occurrence for the entire proteome and dividing by the total number of side chains (total protein length of the entire genome). [This version skips amino acids that are not in the basic 20 (U) and gaps (X). As such the length can be longer than the actual number of side chains counted]

It will parse through all of the GENBANK and FASTA files in a directory and generate two outputs - 
precalculated_avgARSC.txt - this contains the values for each genome and the AVERAGE ARSC values for the genome
individual_protein_ARSC.txt - this contains the counts and lengths of each protein for each genome [LARGE FILE]

CAUTION: Calculating these values for all of the COMPLETE Bacteria genomes (~3,000 genomes, ~7,000 scaffolds, ~8,000,000 proteins) consumed ~10GB RAM and ~72 hrs of computation time.

As such, a precalculated_avgARSC.txt and individual_protein_ARSC.txt has been provided for the COMPLETE Bacteria. The script will process new files (ex. COMPLETE Archaea genomes) and not spend time or memory computing previous genomes.

1. Grzymski and Dussaq. The significance of nitrogen cost minimization in proteomes of marine microorganisms. ISME (2012)
2. Baudouin-Cornu, et al. Molecular evolution of protein atomic composition. Science (2001)

###Dependencies

* [Biopython](http://biopython.org/wiki/Download)

* [pylab](http://www.scipy.org/install.html) through the SciPy package

###Usage
```
python multiGenomeARSC_V2.py
```

##Visualizing the Results
visualizing_genomeavg.py will display a scatterplot of Avg ARSC plotted against %GC and corresponding linear regression for each element

The input is the precalculated_avgARSC.txt file and the output is an HTML file with several functions, such as zooming and panning AND most importantly the ability to hover over single point and see the %GC values, Avg. ARSC values, the species name, and the NCBI Accession number.

###Dependencies

* [Bokeh](http://bokeh.pydata.org/docs/installation.html)

* [Pandas, numpy, and pylab](http://www.scipy.org/install.html) through the SciPy package

* [sklearn](http://scikit-learn.org/stable/install.html)

* [patsy](https://patsy.readthedocs.org/en/latest/overview.html#installation)

###Usage
```
python visualizing_genomeavg.py <FILE NAME | precalculated_avgARSC.txt> <Element of choice: nitrogen, carbon, sulfur, oxygen>
```

#Ongoing Project
Things I would like to expand this project to include:
* Is the correlation seen for carbon and oxygen meaningful?
* Archaeal genomes
* Bacterial and Archaeal draft genomes
* The ability to take a single genome and transition to a visualization of the proteome of that genome
* Use machine learning to determine outliers of interest within the proteome
