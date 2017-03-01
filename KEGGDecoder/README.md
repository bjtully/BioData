KEGG-Decoder
================================================================
Designed to parse through a blastKoala or ghostKoala output to determine the completeness of various KEGG pathways.

This script was constructed using the canonical pathways described as part of KEGG Pathway Maps. There is no additional information provided — if you are interested in certain pathway and the genes are listed in KEGG it is possible to add it to file (with some Python scripting)

###Dependencies

* [Pandas] (http://pandas.pydata.org/pandas-docs/stable/install.html)

* [Seaborn] (http://seaborn.pydata.org/installing.html)

* [matplotlib] (http://matplotlib.org/users/installing.html)

##Additional Information
* Details as to which KEGG KO ids and genes are in each described pathway or process can be found in the supporting document, KOALA_definitions.txt

##Procedure
* Process protein sequences through BlastKOALA or GhostKOALA and download the predicted function text file
* Be sure your submitted FASTA file has headers that group genomes together, KEGG-decoder.py groups based on the name provided in FASTA header before the first underscore (_)
```
For example
>NORP9_1
>NORP9_2
>NORP9_3
>NORP10_1
>NORP10_2
>NORP10_3
Produces two rows of output, one for genome NORP9 and one for genome NORP10 in the list and heat map
```

* The KOALA output text file should look like this:
```
NORP9_1	K00370
NORP9_2	K00371
```
* Run the script
```
python KEGG-decoder.py <INPUT KOALA FILE> <OUTPUT LIST>
```
* The OUTPUT LIST generates a text version of the heat map. The first row contains pathway/process names, subsequent rows contain submitted groups/genomes and fractional percentage of pathway/process

* Figure is output as function_heatmap.svg. Each distinct identifier before the underscore in the FASTA file will have a row