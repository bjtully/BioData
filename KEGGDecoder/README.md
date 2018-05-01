KEGG-Decoder
================================================================
## V0.6 ##
V.0.6 Adds Bacterial Secretion Systems as descrived by KEGG covering Type I, II, III, IV, Vabc, VI, Sec-SRP and Twin Arginine Targeting systems

## V0.5 ##
Adds parameters to force labels to be printed on heatmap. Includes functions
for sulfolipid biosynthesis (key gene sqdB) and C-P lyase

## V0.4 ##
Adds sections that more accurately represents anoxygenic photosynthesis - type-II and type-I reaction centers, adds NiFe hydrogenase Hyd-1 hyaABC, corrected typo leading to missed assignment to hydrogen:quinone oxidoreductase

## V0.3 ##
Latest version adds checks for: retinal biosynthesis, sulfite dehydrogenase (quinone), hydrazine dehydrogenase, hydrazine synthase, DMSP/DMS/DMSO cycling, cobalamin biosynthesis, competence-related DNA transport, anaplerotic reactions

### Description ###
Designed to parse through a blastKoala or ghostKoala output to determine the completeness of various KEGG pathways.

This script was constructed using the canonical pathways described as part of KEGG Pathway Maps. There is no additional information provided — if you are interested in certain pathway and the genes are listed in KEGG it is possible to add it to file (with some Python scripting)

### Dependencies ###

* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html)

* [Seaborn](http://seaborn.pydata.org/installing.html)

* [matplotlib](http://matplotlib.org/users/installing.html)

## Additional Information ##
* Details as to which KEGG KO ids and genes are in each described pathway or process can be found in the supporting document, KOALA_definitions.txt

## Procedure ##
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

KEGG-Expander
================================================================
### Description ###
Designed to expand on the output from KEGG-Decoder. Within KEGG there is a lack of information regarding several processes of interest. To overcome these shortcomings, a small targeted HMM database was created (and will be updated) to fill in gaps of information.

HMM models are predominantly from the PFam database, but when necessary are pulled from TIGRfam and SFam.

### Dependencies ###
* [HMMER3](http://www.hmmer.org/)

### Additional Information ###
* Details as to which HMM models and genes are in each described pathway or process can be found in the supporting document, Pfam_definitions.txt
* In version 0.3, KEGG-Expander targets: phototrophy via proteorhodopsin, peptidases, alternative nitrogenases, ammonia transport, DMSP lyase, and DMSP synthase
* Unfortunately, accuracy depends on the model used, using a bit score cutoff of 75 (approximately an E-value <10E-20) does not always capture the best matches. For example the rhodopsin model does not distinguish between proteorhodopsin and other light driven rhodopsins (we use a tree to determine the proteorhodopsins). Or several of the DMSP lyases at low bit scores will match metalloproteases; in this instance the script has been modified to look for a more stringent bit score (>500). Or the TIGRfam models for the Fe-only and Vanadium nitrogenases generally match the same protein. 

## Prodecure ##
* Using a protein FASTA file with the same gene name set-up as described above - GENOMEID_Number - run a search against the custom HMM database
```
hmmsearch --tblout <NAME>_expanderv0.3.tbl -T 75 /path/to/BioData/KEGGDecoder/HMM_Models/expander_dbv0.3.hmm <INPUT FASTA FILE>
```
* The HMM results table is used to construct the heatmap by running KEGG-expander.py
```
python KEGG-expander.py <NAME>_expanderv0.3.tbl <OUTNAME>_hmm.list
```
* The OUTPUT LIST generates a text version of the heat map. The first row contains pathway/process names, subsequent rows contain submitted groups/genomes and fractional percentage of pathway/process

* Figure is output as hmm_heatmap.svg. Each distinct identifier before the underscore in the FASTA file will have a row

Decoder and Expand
================================================================
### Description ###
Combines the KEGG and HMM heatmaps in to a final heat map. 

### Procedure ###
* Run the script Decoder_and_Expand.py
```
python Decode_and_Expand.py <NAME>_ko.list <NAME>_hmm.list
```
* Figure is output as decode-expand_heatmap.py. Each distinct identifier before the underscore in the FASTA file will have a row