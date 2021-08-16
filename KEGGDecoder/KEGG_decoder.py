#!/usr/bin/python

'''
KEGG-decoder.py V.1.3
V.1.3
Added several pathways associated with carotenoid biosynthesis, including
end-products: astaxanthin, nostoxanthin, zeaxanthin diglucoside, & 
myxoxanthophylls. Plus, staphyloaxanthin biosynthesis and the two pathways for
terpenoid building blocks, the mevalonate pathway and the MEP/DOXP pathway.
The pathways were provided by Dr. Tania Kurbessoian
V.1.2.1
Fixed typo in determing reverse TCA cycle as identified by KEGG-Decoder
user Cheng. Added all-trans-8'-apo-beta-carotenal 15,15'-oxygenase 
which will cleave apo-carotenals to generate retinal. Suggested by Eric Webb.
Upstream pathway unknown
V.1.2
Added several new pathways including PET degradation, carbon storage,
related to starch/gylcogen & polyhydroxybutyrate, and posphate storage,
related to the reversible polyphosphate reaction. Part of summer research
with Sheyla Aviles.
V.1.1
Correcting typos identified by Chris Neely. Adding more complete
pathways components for amino acid biosynthesis identified by
Dr. Eric Webb
phenylalanine added K01713 pheC; cyclohexadienyl dehydratase OR K05359 ADT; 
arogenate/prephenate dehydratase OR K04518 pheA2; prephenate dehydratase
tyrosine added K00220 tyrC; cyclohexadieny/prephenate dehydrogenase OR 
K24018; cyclohexadieny/prephenate dehydrogenase OR K15226 tyrAa; arogenate dehydrogenase
V.1.0.10
Add the biosynthesis of the 20 amino acids - represented as the last
step in the pathway
KEGG-decoder.py V.1.0.8.2
V.1.0.8
Several recent updates have improved all three outputs for visualization
expanded further in the ReadMe note. Additionally, a correction to 
determining the completeness of ubiquinol-cytochrome c reductase. Previously,
only checked for the presence of K00411 and K00410. K00410 is a fusion of
K00412 and K00413 only present in a subset of Proteobacteria. Identified
by Grayson Chadwick
V.1.0.5
Added tanglegram correction for minimizing euclidean distance
V.1.0.4
Removed check for small datasets in interactive viz, added autosizing of tanglegram plots,
switch similarity metric in tanglegram to braycurtis
V.1.0.3
Quality of life updates for KEGG-decoder, tangelgram option, and interactive plots
V.1.0.2
Adds Na+-transporting NADH:ubiquinone oxidoreductase and several metal transporters
KEGG-decoder.py V.0.8
V.0.8
Add elements regarding arsenic reduction
V.0.7
Clarifies elements of methane oxidation and adds additional methanol/alcohol dehydrogenase
to KEGG function search. Adds the serine pathway for formaldehyde assimilation
V.0.6.1 Corrects an issue with the Wood-Ljungdhal pathway that used the wrong
carbon-monoxide deydrogenase subunit
V.0.6 Adds Bacterial Secretion Systems as descrived by KEGG covering Type I, II, III, IV, Vabc,
VI, Sec-SRP and Twin Arginine Targeting systems
V.0.5 Adds parameters to force labels to be printed on heatmap. Includes functions
for sulfolipid biosynthesis (key gene sqdB) and C-P lyase
V.0.4 Adds sections that more accurately represents anoxygenic photosynthesis
- type-II and type-I reaction centers, adds NiFe hydrogenase Hyd-1 hyaABC,
corrected typo leading to missed assignment to hydrogen:quinone oxidoreductase
V.0.3. Adds retinal biosynthesis, sulfite dehydrogenase (quinone), 
hydrazine dehydrogenase, hydrazine synthase, DMSP/DMS/DMSO cycling, 
cobalamin biosynthesis, competence-related DNA transport, anaplerotic 
reactions
Usage: python KEGG-decoder.py <KOALA INPUT> <FUNCTION LIST FORMAT>

Designed to parse through a blastKoala or ghostKoala output to determine
the completeness of various KEGG pathways

Dependencies:
Pandas - http://pandas.pydata.org/pandas-docs/stable/install.html
Seaborn - http://seaborn.pydata.org/installing.html
matplotlib - http://matplotlib.org/users/installing.html

For extended information about KEGG assignments, genes and pathways,
please see accompanying document "KOALA_definitions.txt"

'''
def nitrogen(ko_match):
	out_data = {'dissim nitrate reduction': 0, 'nitrite oxidation': 0,
	'DNRA': 0, 'nitrite reduction': 0, 'nitric oxide reduction' : 0,
	'nitrous-oxide reduction': 0, 'nitrogen fixation' : 0,
	'hydroxylamine oxidation' :0, 'ammonia oxidation (amo/pmmo)': 0,
	'hydrazine dehydrogenase': 0, 'hydrazine synthase': 0}
#narGH
	if ('K00370' in ko_match and 'K00371' in ko_match):
		out_data['dissim nitrate reduction'] = 1
#napAB
	if ('K02567' in ko_match and 'K02568' in ko_match):
		out_data['dissim nitrate reduction'] = 1
#nxrAB
	if ('K00370' in ko_match and 'K00371' in ko_match):
		out_data['nitrite oxidation'] = 1
#nirBD
	if ('K00362' in ko_match and 'K00363' in ko_match):
		out_data['DNRA'] = 1
#nrfAH
	if ('K03385' in ko_match and 'K15876' in ko_match):
		out_data['DNRA'] = 1
#nirK
	if ('K00368' in ko_match):
		out_data['nitrite reduction'] = 1
#nirS
	if ('K15864' in ko_match):
		out_data['nitrite reduction'] = 1
#norBC
	if ('K04561' in ko_match and 'K02305' in ko_match):
		out_data['nitric oxide reduction'] = 1
#nosZ
	if ('K00376' in ko_match):
		out_data['nitrous-oxide reduction'] = 1
#nifKDH
#    if ('K02586' in ko_match and 'K02591' in ko_match and 'K02588' in ko_match):
#        out_data['nitrogen fixation'] = 1
	if ('K02586' in ko_match):
		out_data['nitrogen fixation'] += 0.33
	if ('K02591' in ko_match):
		out_data['nitrogen fixation'] += 0.33
	if ('K02588' in ko_match):
		out_data['nitrogen fixation'] += 0.33
#hao
	if ('K10535' in ko_match):
		out_data['hydroxylamine oxidation'] = 1
#amoA
	if ('K10944' in ko_match):
		out_data['ammonia oxidation (amo/pmmo)'] = 0.33
#amoB
		if ('K10945' in ko_match):
			out_data['ammonia oxidation (amo/pmmo)'] += 0.33
#amoC
		if ('K10946' in ko_match):
			out_data['ammonia oxidation (amo/pmmo)'] += 0.33
	if ('K20935' in ko_match):
		out_data['hydrazine dehydrogenase'] = 1
	hydrazine_synth = ['K20932', 'K20933', 'K20934']
	for i in hydrazine_synth:
		if i in ko_match:
			out_data['hydrazine synthase'] += 0.33

	return out_data

def glycolysis(ko_match):
#Check for presence of 9 genes
	total = 0
#phosphoglucomutase, glucose-6-phosphate isomerase, fructose-bisphosphate aldolase
#phosphoglycerate kinase, enolase
	single_ko = ['K01835', 'K01810', 'K01623', 'K00927', 'K01689']
	for i in single_ko:
		if i in ko_match:
			total += 1
#6-phosphofructokinase
	if ('K00850' in ko_match or 'K00895' in ko_match):
		total += 1
#glyceraldehyde 3-phosphate dehydrogenase
	if ('K00134' in ko_match or 'K00150' in ko_match):
		total += 1
#2,3-bisphosphoglycerate-dependent phosphoglycerate mutase 
	if ('K01834' in ko_match or 'K15633' in ko_match):
		total += 1
#pyruvate kinase
	if ('K00873' in ko_match or 'K01006' in ko_match):
		total += 1
	value = float(total)/float(9)
	return {'glycolysis': float("%.2f" % (value))}

def gluconeogenesis(ko_match):
	total = 0
#Requires fructose-1,6-bisphosphatase to continue
	if ('K03841' in ko_match):
		total += 1
#phosphoglucomutase, glucose-6-phosphate isomerase, fructose-bisphosphate aldolase
#phosphoglycerate kinase, enolase
		single_ko = ['K01835', 'K01810', 'K01623', 'K00927', 'K01689']
		for i in single_ko:
			if i in ko_match:
				total += 1
#glyceraldehyde 3-phosphate dehydrogenase
		if ('K00134' in ko_match or 'K00150' in ko_match):
			total += 1
#2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
		if ('K01834' in ko_match or 'K15633' in ko_match):
			total += 1
#pyruvate kinase
		if ('K00873' in ko_match or 'K01006' in ko_match):
			total += 1
	value = float(total)/float(9)
	return {'gluconeogenesis': float("%.2f" % (value))}

def tca_cycle(ko_match):
	total = 0
#aconitate hydratase
	if ('K01681' in ko_match or 'K01682' in ko_match):
		total += 1
#isocitrate dehydrogenase
	if ('K00031' in ko_match) or ('K00030' in ko_match) or ('K17753' in ko_match):
		total += 1
#2-oxoglutarate/2-oxoacid ferredoxin oxidoreductase
	if ('K00174' in ko_match and 'K00175' in ko_match):
		total += 1
#succinyl-CoA synthetase
	if (('K01899' in ko_match and 'K01900' in ko_match) or
		('K01902' in ko_match and 'K01903' in ko_match) or
		('K18118' in ko_match)):
		total += 1
#fumarate reductase 
	if (('K00244' in ko_match and 'K00245' in ko_match and 'K00246' in ko_match and 'K00247' in ko_match)
		or
		('K00239' in ko_match and 'K00240' in ko_match and 'K00241' in ko_match and 'K00242' in ko_match)
		or
		('K00234' in ko_match and 'K00235' in ko_match and 'K00236' in ko_match and 'K00237' in ko_match)):
		total += 1
#fumurate hydratase
	if (('K01677' in ko_match and 'K01678' in ko_match and 'K01679' in ko_match) or
		('K01676' in ko_match)):
		total += 1
#malate dehydrogenase
	if (('K00116' in ko_match) or
		('K00025' in ko_match) or
		('K00026' in ko_match) or
		('K00024' in ko_match)):
		total += 1
#citrate synthase
	if ('K01647' in ko_match):
		total += 1
	value = float(total)/float(8)
	return {'TCA Cycle': float("%.2f" % (value))}

def cbb_cycle(ko_match):
	total = 0
	var_cnt = 4
	out_data = {'RuBisCo' : 0, 'CBB Cycle': 0}
#RuBisCO - Only large subunit Type 1 and 2
	if ('K01601' in ko_match):
		out_data['RuBisCo'] = 1
		total += 1
#phosphoglycerate kinase
		if ('K00927' in ko_match):
			total += 1
#glyceraldehyde 3-phosphate dehydrogenase 
		if ('K00134' in ko_match) or ('K05298' in ko_match) or ('K00150' in ko_match):
			total += 1
#phosphoribulokinase            
		if ('K00855' in ko_match):
			total += 1
#Ribulose regeneration
#ribulose-phosphate 3-epimerase AND xylulose-5-phosphate/fructose-6-phosphate phosphoketolase
		if ('K01783' in ko_match and 'K01621' in ko_match):
			total += 2
			var_cnt += 2
#transketolase AND ribulose-phosphate 3-epimerase
		if ('K00615' in ko_match and 'K01783' in ko_match):
			total += 2
			var_cnt += 2
#transketolase AND ribose 5-phosphate isomerase
		if ('K00615' in ko_match and 'K01807' in ko_match):
			total += 2
			var_cnt += 2
#fructose-bisphosphate aldolase AND transketolase AND fructose-1,6-bisphosphatase
		if (('K01623' in ko_match or 'K01624' in ko_match or 'K11645' in ko_match) and
			('K00615' in ko_match) and
			('K11532' in ko_match or 'K03841' in ko_match or 'K02446' in ko_match)):
			total += 3
			var_cnt += 3
	value = float(total)/float(var_cnt)
	out_data['CBB Cycle'] = float("%.2f" % (value))
	return out_data

def reverse_tca(ko_match):
	out_data = {'rTCA Cycle' : 0}
#ATP-citrate lyase
	if ('K15230' in ko_match and 'K15231' in ko_match):
		out_data['rTCA Cycle'] = 1
#citryl-CoA synthetase AND citryl-CoA lyase
	if ('K15232' in ko_match and 'K15233' in ko_match and 'K15234' in ko_match):
		out_data['rTCA Cycle'] = 1
	return out_data

def wood_ljungdahl(ko_match):
	total = 0
	CO_methyl_present = 0
#Carbon fixing branch
#acetyl-CoA decarbonylase/synthase complex subunit alpha OR 
#CO-methylating acetyl-CoA synthase
	if ('K00192' in ko_match) or ('K14138' in ko_match):
		total += 1
		CO_methyl_present = 1
#catalytic subunits only of CO dehydrogenase
#anaerobic carbon-monoxide dehydrogenase OR aerobic carbon-monoxide dehydrogenase large subunit
	if ('K00198' in ko_match) or ('K03520' in ko_match):
		total+= 1
	if CO_methyl_present == 1:
#Methyl branch
#formate dehydrogenase
		if ('K05299' in ko_match and 'K15022' in ko_match):
			total+= 1
#formate--tetrahydrofolate ligase
		if ('K01938' in ko_match):
			total+= 1
#methylenetetrahydrofolate dehydrogenase (NADP+) / methenyltetrahydrofolate cyclohydrolase
		if ('K01491' in ko_match):
			total+= 1
#methylenetetrahydrofolate reductase (NADPH)
		if ('K00297' in ko_match):
			total+= 1
	value = float(total)/float(6)
	return {'Wood-Ljungdahl' : float("%.2f" % (value))}

def three_prop(ko_match):
	total = 0
#pyruvate ferredoxin oxidoreductase alpha and beta subunits
	if ('K00169' in ko_match and 'K00170' in ko_match):
		total +=1
#pyruvate dikinase
	if ('K01006' in ko_match or 'K01007' in ko_match):
		total +=1
#phosphoenolpyruvate carboxylase
	if ('K01595' in ko_match):
		total +=1
#malate dehydrogenase
	if ('K00024' in ko_match):
		total +=1
#succinyl-CoA:(S)-malate CoA-transferase
	if ('K14471' in ko_match and 'K14472' in ko_match):
		total +=1
#malyl-CoA/(S)-citramalyl-CoA lyase
	if ('K08691' in ko_match):
		total +=1
#acetyl-CoA carboxylase, biotin carboxylase
	if ('K02160' in ko_match and 'K01961' in ko_match and
		'K01962' in ko_match and 'K01963' in ko_match):
		total +=1
#malonyl-CoA reductase / 3-hydroxypropionate dehydrogenase (NADP+) 
	if ('K14468' in ko_match and 'K15017' in ko_match):
		total +=1
#3-hydroxypropionate dehydrogenase (NADP+)
	if ('K15039' in ko_match):
		total +=1
#acrylyl-CoA reductase (NADPH) / 3-hydroxypropionyl-CoA dehydratase / 3-hydroxypropionyl-CoA synthetase
	if ('K14469' in ko_match and 'K15018' in ko_match):
		total +=1
#3-hydroxypropionyl-coenzyme A dehydratase
	if ('K15019' in ko_match):
		total +=1
#acryloyl-coenzyme A reductase
	if ('K15020' in ko_match):
		total +=1
#malyl-CoA/(S)-citramalyl-CoA lyase
	if ('K08691' in ko_match):
		total +=1
#2-methylfumaryl-CoA hydratase
	if ('K14449' in ko_match):
		total +=1
#2-methylfumaryl-CoA isomerase 
	if ('K14470' in ko_match):
		total +=1
#3-methylfumaryl-CoA hydratase
	if ('K09709' in ko_match):
		total +=1
#malyl-CoA/(S)-citramalyl-CoA lyase
	if ('K08691' in ko_match):
		total +=1
	value = float(total)/float(17)
	return {'3-Hydroxypropionate Bicycle' : float("%.2f" % (value))}

def four_hydrox(ko_match):
#Based on the reference present in Thaumarchaea -- pathway is not complete
	total = 0
#acetyl-CoA carboxylase, biotin carboxylase
	if ('K02160' in ko_match and 'K01961' in ko_match and
		'K01962' in ko_match and 'K01963' in ko_match):
		total +=1
#malonic semialdehyde reductase
	if ('K18602' in ko_match):
		total +=1
#3-hydroxypropionyl-CoA synthetase
	if ('K18594' in ko_match):
		total +=1
#acrylyl-CoA reductase (NADPH) / 3-hydroxypropionyl-CoA dehydratase / 3-hydroxypropionyl-CoA synthetase
	if ('K14469' in ko_match and 'K15019' in ko_match):
		total +=1
#methylmalonyl-CoA/ethylmalonyl-CoA epimerase
	if ('K05606' in ko_match):
		total +=1
#methylmalonyl-CoA mutase
	if ('K01847' in ko_match and 'K01848' in ko_match and
		'K01849' in ko_match):
		total +=1
#4-hydroxybutyryl-CoA synthetase (ADP-forming)
	if ('K18593' in ko_match):
		total +=1
#4-hydroxybutyryl-CoA dehydratase / vinylacetyl-CoA-Delta-isomerase
	if ('K14534' in ko_match):
		total +=1
#enoyl-CoA hydratase / 3-hydroxyacyl-CoA dehydrogenase
	if ('K15016' in ko_match):
		total +=1
#acetyl-CoA C-acetyltransferase
	if ('K00626' in ko_match):
		total +=1
	value = float(total)/float(10)
	return {'4-Hydroxybutyrate/3-hydroxypropionate' : float("%.2f" % (value))}

def c_degradation(ko_match):
	out_data = {'beta-glucosidase' : 0, 'cellulase':0, 'chitinase':0,
				'bifunctional chitinase/lysozyme':0,
				'basic endochitinase B':0, 'diacetylchitobiose deacetylase':0,
				'beta-N-acetylhexosaminidase':0, 'pectinesterase':0,
				'exo-poly-alpha-galacturonosidase':0, 'oligogalacturonide lyase':0,
				'exopolygalacturonase':0, 'D-galacturonate isomerase':0,
				'D-galacturonate epimerase':0, 'alpha-amylase': 0, 'glucoamylase':0,
				'pullulanase':0}
	if ('K05350' in ko_match or 'K05349' in ko_match):
		out_data['beta-glucosidase'] += 1
	if ('K01225' in ko_match or 'K19668' in ko_match):
		out_data['cellulase'] += 1
	if ('K01183' in ko_match):
		out_data['chitinase'] += 1
	if ('K13381' in ko_match):
		out_data['bifunctional chitinase/lysozyme'] += 1
	if ('K20547' in ko_match):
		out_data['basic endochitinase B'] += 1
	if ('K03478' in ko_match or 'K18454' in ko_match):
		out_data['diacetylchitobiose deacetylase'] += 1
	if ('K01207' in ko_match):
		out_data['beta-N-acetylhexosaminidase'] += 1
	if ('K01730' in ko_match):
		out_data['oligogalacturonide lyase'] += 1
	if ('K01184' in ko_match):
		out_data['exopolygalacturonase'] += 1
	if ('K01812' in ko_match):
		out_data['D-galacturonate isomerase'] += 1
	if ('K08679' in ko_match):
		out_data['D-galacturonate epimerase'] += 1
	if ('K01176' in ko_match):
		out_data['alpha-amylase'] += 1
	if ('K01178' in ko_match):
		out_data['glucoamylase'] += 1
	if ('K01200' in ko_match):
		out_data['pullulanase'] += 1
	return out_data

def chemotaxis(ko_match):
#Che family of proteins
	total = 0
	single_ko = ['K13924', 'K00575',
				'K03413', 'K03412', 'K03406', 'K03407',
				'K03415', 'K03408']
	for i in single_ko:
		if i in ko_match:
			total += 1
	value = float(total)/float(len(single_ko))
	return {'Chemotaxis': float("%.2f" % (value))}

def flagellum(ko_match):
#Components of the flagellum biosynthesis group
	total = 0
	single_ko = ['K02409', 'K02401', 'K02394', 'K02397',
				'K02396', 'K02391', 'K02390', 'K02393',
				'K02392', 'K02386', 'K02557', 'K02556',
				'K02400', 'K02418', 'K02389', 'K02412',
				'K02387', 'K02410', 'K02411', 'K02416',
				'K02417', 'K02407', 'K02406']
	for i in single_ko:
		if i in ko_match:
			total += 1
	value = float(total)/float(len(single_ko))
	return {'Flagellum': float("%.2f" % (value))}

def sulfur(ko_match):
	out_data = {'sulfur assimilation':0, 'dissimilatory sulfate < > APS':0,
		'dissimilatory sulfite < > APS':0, 'dissimilatory sulfite < > sulfide':0,
		'thiosulfate oxidation':0, 'alt thiosulfate oxidation doxAD':0,
		'alt thiosulfate oxidation tsdA':0, 'thiosulfate disproportionation':0, 'sulfur reductase sreABC':0,
		'thiosulfate/polysulfide reductase':0, 'sulfhydrogenase':0,
		'sulfur disproportionation':0, 'sulfur dioxygenase':0, 'sulfite dehydrogenase':0,
		'sulfide oxidation':0, 'sulfite dehydrogenase (quinone)':0,
		'DMSP demethylation': 0, 'DMS dehydrogenase': 0, 'DMSO reductase': 0}
#sir; sulfite reductase (ferredoxin) [EC:1.8.7.1] OR
#cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2] + cysI; sulfite reductase (NADPH) hemoprotein beta-component [EC:1.8.1.2]
	if ('K00392' in ko_match) or ('K00380' in ko_match and 'K00381' in ko_match):
		out_data['sulfur assimilation'] = 1
#sat; sulfate adenylyltransferase
	if ('K00958' in ko_match):
		out_data['dissimilatory sulfate < > APS'] = 1
#aprB; adenylylsulfate reductase, subunit B [EC:1.8.99.2] + aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]
	if ('K00395' in ko_match and 'K00394' in ko_match):
		out_data['dissimilatory sulfite < > APS'] = 1
#dsrA; sulfite reductase, dissimilatory-type alpha subunit [EC:1.8.99.3] + dsrB; sulfite reductase, dissimilatory-type beta subunit [EC:1.8.99.3]
	if ('K11180' in ko_match and 'K11181' in ko_match):
		out_data['dissimilatory sulfite < > sulfide'] = 1
#soxABCXYZ
	if ('K17222' in ko_match):
		out_data['thiosulfate oxidation'] += .16
	if ('K17224' in ko_match):
		out_data['thiosulfate oxidation'] += .16
	if ('K17225' in ko_match):
		out_data['thiosulfate oxidation'] += .16
	if ('K17223' in ko_match):
		out_data['thiosulfate oxidation'] += .16
	if ('K17226' in ko_match):
		out_data['thiosulfate oxidation'] += .16
	if ('K17227' in ko_match):
		out_data['thiosulfate oxidation'] += .16
#doxAD thiosulfate dehydrogenase [quinone]
	if ('K16936' in ko_match and 'K16937' in ko_match):
		out_data['alt thiosulfate oxidation doxAD'] = 1
#tsdA thiosulfate dehydrogenase [EC:1.8.2.2]
	if ('K19713' in ko_match):
		out_data['alt thiosulfate oxidation tsdA'] = 1
#sulfur reductase 
	if ('K17219' in ko_match):
		out_data['sulfur reductase sreABC'] += .33
	if ('K17220' in ko_match):
		out_data['sulfur reductase sreABC'] += .33
	if ('K17221' in ko_match):
		out_data['sulfur reductase sreABC'] += .33
#thiosulfate reductase / polysulfide reductase psrABC/phsABC
	if ('K08352' in ko_match):
		out_data['thiosulfate/polysulfide reductase'] += .33
	if ('K08353' in ko_match):
		out_data['thiosulfate/polysulfide reductase'] += .33
	if ('K08354' in ko_match):
		out_data['thiosulfate/polysulfide reductase'] += .33
#sulfhydrogenase hydABGD
	if ('K17993' in ko_match):
		out_data['sulfhydrogenase'] += .25
	if ('K17996' in ko_match):
		out_data['sulfhydrogenase'] += .25
	if ('K17995' in ko_match):
		out_data['sulfhydrogenase'] += .25
	if ('K17994' in ko_match):
		out_data['sulfhydrogenase'] += .25
#sor; sulfur oxygenase/reductase
	if ('K16952' in ko_match):
		out_data['sulfur disproportionation'] += 1
#sdo; sulfur dioxygenase
	if ('K17725' in ko_match):
		out_data['sulfur dioxygenase'] += 1
#sorB; sulfite dehydrogenase
	if ('K05301' in ko_match):
		out_data['sulfite dehydrogenase'] += 1
#sqr; sulfide:quinone oxidoreductase OR fccB; sulfide dehydrogenase [flavocytochrome c]
	if ('K17218' in ko_match or 'K17229' in ko_match):
		out_data['sulfide oxidation'] = 1
	value = out_data['thiosulfate oxidation']
	out_data['thiosulfate oxidation'] = float("%.2f" % (value))
#soeABC; sulfite dehydrogenase (quinone) 
	soeABC = ['K21307', 'K21308', 'K21309']
	for i in soeABC:
		if i in ko_match:
			out_data['sulfite dehydrogenase (quinone)'] += 0.33
#DMSP lyase
#    if ('K16953' in ko_match):
#        out_data['DMSP lyase, dddL'] = 1
#dmdA; dimethylsulfoniopropionate demethylase
	if ('K17486' in ko_match):
		out_data['DMSP demethylation'] = 1
#ddhABC; dimethylsulfide dehydrogenase
	dms_dh = ['K16964', 'K16965', 'K16966']
	for i in dms_dh:
		if i in ko_match:
			out_data['DMS dehydrogenase'] += 0.33
#dmsABC; anaerobic dimethyl sulfoxide reductase
	dmso_red = ['K07306', 'K07307', 'K07308']
	for i in dmso_red:
		if i in ko_match:
			out_data['DMSO reductase'] += 0.33
	return out_data

def methanogenesis(ko_match):
	out_data = {'Methanogenesis via methanol':0, 'Methanogenesis via dimethylamine':0,
		'Methanogenesis via dimethylsulfide, methanethiol, methylpropanoate':0,
		'Methanogenesis via methylamine':0, 'Methanogenesis via trimethylamine':0,
		'Methanogenesis via acetate':0, 'Methanogenesis via CO2':0,
		'Coenzyme M reduction to methane':0, 'Coenzyme B/Coenzyme M regeneration':0,
		'dimethylamine/trimethylamine dehydrogenase':0}
#dmd-tmd; dimethylamine/trimethylamine dehydrogenase
	if ('K00317' in ko_match):
		out_data['dimethylamine/trimethylamine dehydrogenase'] = 1
#mtaA; [methyl-Co(III) methanol-specific corrinoid protein]:coenzyme M methyltransferase
#mtaB; methanol---5-hydroxybenzumidazolylcobamide Co-methyltransferase
#mtaC; methanol corrinoid protein
	methanol_kos = ['K14080', 'K04480', 'K14081']
	for i in methanol_kos:
		if i in ko_match:
			out_data['Methanogenesis via methanol'] += .33
#mtbA; [methyl-Co(III) methylamine-specific corrinoid protein]:coenzyme M methyltransferase
#mtbB; dimethylamine---corrinoid protein Co-methyltransferase
	dimethylamine_kos = ['K14082', 'K16178']
	for i in dimethylamine_kos:
		if i in ko_match:
			out_data['Methanogenesis via dimethylamine'] += .50
#mtsA; methylthiol:coenzyme M methyltransferase
#mtsB; methylated-thiol--corrinoid protein
	dimethylsulfide_kos = ['K16954', 'K16955']
	for i in dimethylsulfide_kos:
		if i in ko_match:
			out_data['Methanogenesis via dimethylsulfide, methanethiol, methylpropanoate'] += .50
#mtmB; monomethylamine methyltransferase
	if ('K16178' in ko_match):
		out_data['Methanogenesis via methylamine'] = 1
#mttB; trimethylamine methyltransferase
	if ('K14083' in ko_match):
		out_data['Methanogenesis via trimethylamine'] = 1
#acetyl-CoA decarbonylase/synthase complex
	acetate_kos = ['K00193', 'K00194', 'K00197']
	for i in acetate_kos:
		if i in ko_match:
			out_data['Methanogenesis via acetate'] += .33
#formylmethanofuran dehydrogenase
#ftr; formylmethanofuran--tetrahydromethanopterin N-formyltransferase
#mch; methenyltetrahydromethanopterin cyclohydrolase
#hmd; 5,10-methenyltetrahydromethanopterin hydrogenase
#mer; 5,10-methylenetetrahydromethanopterin reductase
#mtrABCDEFGH; tetrahydromethanopterin S-methyltransferase
	co2_kos = ['K00200', 'K00201', 'K00202', 'K00203', 'K00205', 'K11261',
		'K00672', 'K01499', 'K13942', 'K00320', 'K00577', 'K00578', 'K00579',
		'K00580', 'K00581', 'K00582', 'K00583', 'K00584']
	for i in co2_kos:
		if i in ko_match:
			out_data['Methanogenesis via CO2'] += .05
#mcrABCD; methyl-coenzyme M reductase
	coenzymeM_kos = ['K00399', 'K00401', 'K00402']
	for i in coenzymeM_kos:
		if i in ko_match:
			out_data['Coenzyme M reduction to methane'] += .33
#hdrABCDE; CoB-CoM heterodisulfide reductase
	regeneration_kos = ['K03388', 'K03389', 'K03390', 'K08264', 'K08265']
	for i in regeneration_kos:
		if i in ko_match:
			out_data['Coenzyme B/Coenzyme M regeneration'] += .20
	value = out_data['Coenzyme B/Coenzyme M regeneration']
	out_data['Coenzyme B/Coenzyme M regeneration'] = float("%.2f" % (value))
	return out_data

def methane_ox(ko_match):
	out_data = {'Soluble methane monooxygenase':0,'methanol dehydrogenase':0,
	'alcohol oxidase':0}
#mmoXYZC; soluble methane monooxygenase
	single_ko = ['K16157', 'K16158', 'K16159', 'K16161']
	for i in single_ko:
		if i in ko_match:
			out_data['Soluble methane monooxygenase'] += .25
	methanol_dh = ['K14028', 'K14029']
	for i in methanol_dh:
		if i in ko_match:
			out_data['methanol dehydrogenase'] += .5
	if ('K17066' in ko_match):
		out_data['alcohol oxidase'] = 1
	return out_data

def hydrogen(ko_match):
	out_data = {'NiFe hydrogenase':0, 'membrane-bound hydrogenase':0,
	'ferredoxin hydrogenase':0, 'hydrogen:quinone oxidoreductase':0,
	'NAD-reducing hydrogenase':0, 'NADP-reducing hydrogenase':0,
	'NiFe hydrogenase Hyd-1':0}
#hydB2,hydA2; NiFe hydrogenase
	if ('K00437' in ko_match and 'K18008' in ko_match):
		out_data['NiFe hydrogenase'] = 1
#mbhLJK; membrane-bound hydrogenase
	if ('K18016' in ko_match and 'K18017' in ko_match
		and 'K18023' in ko_match):
		out_data['membrane-bound hydrogenase'] = 1
#hupSL; ferrodoxin hydrogenase
	if ('K00533' in ko_match and 'K00534' in ko_match):
		out_data['ferredoxin hydrogenase'] = 1
#hydA3,hydB3; hydrogen:quinone oxidoreductase
	if ('K05922' in ko_match and 'K05927' in ko_match):
		out_data['hydrogen:quinone oxidoreductase'] = 1
#hoxHFUY; NAD-reducing hydrogenase
	nad_ko = ['K00436' , 'K18005' , 'K18006' , 'K18007']
	for i in nad_ko:
		if i in ko_match:
			out_data['NAD-reducing hydrogenase'] += .25
#hndABCD; NADP-reducing hydrogenase
	nadp_ko = ['K17992', 'K18330', 'K18331', 'K18332']
	for i in nadp_ko:
		if i in ko_match:
			out_data['NADP-reducing hydrogenase'] += .25
#hyaABC; NiFe hydrogenase Hyd-1
	hyd_ko = ['K06282', 'K06281', 'K03620']
	for i in hyd_ko:
		if i in ko_match:
			out_data['NiFe hydrogenase Hyd-1'] += 0.33
	return out_data

def transporters(ko_match):
	out_data = {'transporter: phosphate':0, 'transporter: phosphonate':0,
	'transporter: thiamin':0, 'transporter: vitamin B12':0,
	'transporter: urea':0}
#pstABCS; phosphate
	phosphate_ko = ['K02040', 'K02037', 'K02038', 'K02036']
	for i in phosphate_ko:
		if i in ko_match:
			out_data['transporter: phosphate'] += .25
#phnDEC; phosphonate
	phosphonate_ko = ['K02044', 'K02042', 'K02041']
	for i in phosphonate_ko:
		if i in ko_match:
			out_data['transporter: phosphonate'] += .33
#tbpA,thiPQ; thiamin
	thiamin_ko = ['K02064', 'K02063', 'K02062']
	for i in thiamin_ko:
		if i in ko_match:
			out_data['transporter: thiamin'] += .33
#btuFCD; vitamin B12
	b12_ko = ['K06858', 'K06073', 'K06074']
	for i in b12_ko:
		if i in ko_match:
			out_data['transporter: vitamin B12'] += .33
#urtABCED; urea
	urea_ko = ['K11959', 'K11960', 'K11961', 'K11962', 'K11963']
	for i in urea_ko:
		if i in ko_match:
			out_data['transporter: urea'] += .2
	value = out_data['transporter: urea']
	out_data['transporter: urea'] = float("%.2f" % (value))
	return out_data

def riboflavin(ko_match):
	total= 0
#ribB; 3,4-dihydroxy 2-butanone 4-phosphate synthase OR 
#ribAB; 3,4-dihydroxy 2-butanone 4-phosphate synthase / GTP cyclohydrolase II
	if ('K02858' in ko_match or 'K14652' in ko_match):
		total += 1
#ribD2; 5-amino-6-(5-phosphoribosylamino)uracil reductase OR 
#ribD; diaminohydroxyphosphoribosylaminopyrimidine deaminase / 5-amino-6-(5-phosphoribosylamino)uracil reductase
	if ('K00082' in ko_match or 'K11752' in ko_match):
		total += 1
#ribH; 6,7-dimethyl-8-ribityllumazine synthase
	if ('K00794' in ko_match):
		total += 1
#ribE; riboflavin synthase
	if ('K00793' in ko_match):
		total += 1
#RFK; riboflavin kinase OR FHY; riboflavin kinase / FMN hydrolase OR 
#ribF; riboflavin kinase / FMN adenylyltransferase
#    if ('K00861' in ko_match or 'K20884' in ko_match or 'K11753' in ko_match):
#        total += 1
#FAD synthetase
#    if ('K14656' in ko_match or 'K00953' in ko_match):
#        total += 1
	value = float(total)/float(4)
	return {'riboflavin biosynthesis': float("%.2f" % (value))}

def thiamin(ko_match):
	total = 0
#thiF; sulfur carrier protein ThiS adenylyltransferase
	if ('K03148' in ko_match):
		total += 1
#iscS; cysteine desulfurase
	if ('K04487' in ko_match):
		total += 1
#thiH; 2-iminoacetate synthase OR K03153 thiO; glycine oxidase
	if ('K03150' in ko_match or 'K03153' in ko_match):
		total += 1
#thiI; thiamine biosynthesis protein ThiI
	if ('K03151' in ko_match):
		total += 1
#dxs; 1-deoxy-D-xylulose-5-phosphate synthase
	if ('K01662' in ko_match):
		total += 1
#thiG; thiazole synthase
	if ('K03149' in ko_match):
		total += 1
#tenI; thiazole tautomerase OR THI4; thiamine thiazole synthase
	if ('K10810' in ko_match or 'K03146' in ko_match):
		total += 1
#THI5; pyrimidine precursor biosynthesis enzyme OR K03147 thiC; phosphomethylpyrimidine synthase OR 
#THI20; hydroxymethylpyrimidine/phosphomethylpyrimidine kinase / thiaminase OR 
#thiD; hydroxymethylpyrimidine/phosphomethylpyrimidine kinase OR 
#thiDE; hydroxymethylpyrimidine kinase / phosphomethylpyrimidine kinase / thiamine-phosphate diphosphorylase
	if ('K18278' in ko_match or 'K03147' in ko_match or
		'K00877' in ko_match or 'K00941' in ko_match):
		total += 1
#THI20; hydroxymethylpyrimidine/phosphomethylpyrimidine kinase / thiaminase OR 
#thiD; hydroxymethylpyrimidine/phosphomethylpyrimidine kinase
	if ('K00877' in ko_match or 'K00941' in ko_match):
		total += 1
#thiE; thiamine-phosphate pyrophosphorylase OR 
#thiDE; hydroxymethylpyrimidine kinase / phosphomethylpyrimidine kinase / thiamine-phosphate diphosphorylase OR 
#THI6;thiamine-phosphate diphosphorylase / hydroxyethylthiazole kinase
	if ('K00788' in ko_match or 'K14153' in ko_match or
		'K14154' in ko_match):
		total += 1
#thiL; thiamine-monophosphate kinase
	if ('K00946' in ko_match):
		total += 1
	value = float(total)/float(11)
	return {'thiamin biosynthesis': float("%.2f" % (value))}

def cobalamin(ko_match):
	total = 0
#pduO; cob(I)alamin adenosyltransferase
#cobA; cob(I)alamin adenosyltransferase
	if ('K00798' in ko_match or 'K19221' in ko_match):
		total += 1
#cobQ; adenosylcobyric acid synthase
	if ('K02232' in ko_match):
		total += 1
#cobC; cobalamin biosynthetic protein CobC
	if ('K02225' in ko_match):
		total += 1
#cobD; adenosylcobinamide-phosphate synthase
	if ('K02227' in ko_match):
		total += 1
#cobU; adenosylcobinamide kinase / adenosylcobinamide-phosphate guanylyltransferase
	if ('K02231' in ko_match):
		total += 1
#cobY; adenosylcobinamide-phosphate guanylyltransferase
	if ('K19712' in ko_match and 'K02231' not in ko_match):
		total += 1
#cobV; adenosylcobinamide-GDP ribazoletransferase
	if ('K02233' in ko_match):
		total += 1
#cobC; alpha-ribazole phosphatase
	if ('K02226' in ko_match):
		total += 1
#cobT; nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase
	if ('K00768' in ko_match):
		total += 1
	value = float(total)/float(8)
	return {'cobalamin biosynthesis': float("%.2f" % (value))}



def oxidative_phoshorylation(ko_match):
	out_data ={'F-type ATPase':0, 'V-type ATPase':0, 'NADH-quinone oxidoreductase':0,
	'NAD(P)H-quinone oxidoreductase':0, 'Cytochrome c oxidase, cbb3-type':0,
	'Cytochrome bd complex':0, 'Cytochrome o ubiquinol oxidase':0,
	'Cytochrome c oxidase':0, 'Cytochrome aa3-600 menaquinol oxidase':0,
	'Ubiquinol-cytochrome c reductase':0, 'Na-NADH-ubiquinone oxidoreductase': 0}
#atpFBCHGDAE
	ftype_ko = ['K02111', 'K02112', 'K02115', 'K02113',
	'K02114', 'K02108', 'K02109', 'K02110']
	for i in ftype_ko:
		if i in ko_match:
			out_data['F-type ATPase'] += 0.125
#ntpABCDEFIK,ahaH
	vtype_ko = ['K02117', 'K02118', 'K02119', 'K02120',
	'K02121', 'K02122', 'K02107', 'K02123', 'K02124']
	for i in vtype_ko:
		if i in ko_match:
			out_data['V-type ATPase'] += 0.11
#nuoABCDEFGHIJKLMN
	nuo_ko = ['K00330', 'K00331', 'K00332', 'K00333',
	'K00334', 'K00335', 'K00336', 'K00337',
	'K00338', 'K00339', 'K00340', 'K00341',
	'K00342', 'K00343']
	for i in nuo_ko:
		if i in ko_match:
			out_data['NADH-quinone oxidoreductase'] += 0.07
	value = out_data['NADH-quinone oxidoreductase']
	out_data['NADH-quinone oxidoreductase'] = float("%.2f" % (value))
#ndcABCDEFGHIJKLMN
	ndc_ko = ['K05574', 'K05582', 'K05581', 'K05579',
	'K05572', 'K05580', 'K05578', 'K05576',
	'K05577', 'K05575', 'K05573', 'K05583',
	'K05584', 'K05585']
	for i in ndc_ko:
		if i in ko_match:
			out_data['NAD(P)H-quinone oxidoreductase'] += 0.07
#ccoPQNO
	cbb3_ko = ['K00404', 'K00405', 'K00407', 'K00406']
	for i in cbb3_ko:
		if i in ko_match:
			out_data['Cytochrome c oxidase, cbb3-type'] += 0.25
#cydAB
	bd_ko = ['K00425', 'K00426']
	for i in bd_ko:
		if i in ko_match:
			out_data['Cytochrome bd complex'] += 0.5
#cyoABCD
	o_ko = ['K02300', 'K02299', 'K02298', 'K02297']
	for i in o_ko:
		if i in ko_match:
			out_data['Cytochrome o ubiquinol oxidase'] += 0.25
#coxABCD
	cytc_ko = ['K02277', 'K02276', 'K02274', 'K02275']
	for i in cytc_ko:
		if i in ko_match:
			out_data['Cytochrome c oxidase'] += 0.25
#qoxABCD
	aa3_ko = ['K02829', 'K02828', 'K02827', 'K02826']
	for i in aa3_ko:
		if i in ko_match:
			out_data['Cytochrome aa3-600 menaquinol oxidase'] += 0.25
#petA,fbcH; ubiquinol-cytochrome c reductase 
#petA,petB,petC; ubiquinol-cytochrome c reductase
	if ('K00411' in ko_match) and ('K00410' in ko_match):
		out_data['Ubiquinol-cytochrome c reductase'] = 1
	else:
		ubiquinol_ko = ['K00411', 'K00412', 'K00413']
		for i in ubiquinol_ko:
			if i in ko_match:
				out_data['Ubiquinol-cytochrome c reductase'] += 0.33
	
#nqrABCDEF; Na+-transporting NADH:ubiquinone oxidoreductase
	na_ubiquinone_ko = ['K00346', 'K00347', 'K00348', 'K00349',
	'K00350', 'K00351']
	for i in na_ubiquinone_ko:
		if i in ko_match:
			out_data['Na-NADH-ubiquinone oxidoreductase'] += 0.167
	value = out_data['Na-NADH-ubiquinone oxidoreductase']
	out_data['Na-NADH-ubiquinone oxidoreductase'] = float("%.2f" % (value))
	
	return out_data

def photosynthesis(ko_match):
	out_data = {'Photosystem II':0, 'Photosystem I':0, 'Cytochrome b6/f complex':0,
	"anoxygenic type-II reaction center":0, "anoxygenic type-I reaction center":0,
	'Retinal biosynthesis':0, 'Retinal from apo-carotenals': 0}
	psII = ['K02703', 'K02706', 'K02705', 'K02704', 'K02707', 'K02708']
#Photosystem II core complex
	for i in psII:
		if i in ko_match:
			out_data['Photosystem II'] += 0.167
	psI = ['K02689', 'K02690', 'K02691', 'K02692', 'K02693', 'K02694', 'K02696',
			'K02697', 'K02698', 'K02699', 'K02700', 'K08905', 'K02695', 'K02701',
			'K14332', 'K02702']
#Photosystem I
	for i in psI:
		if i in ko_match:
			out_data['Photosystem I'] += 0.0625
	cyt_b6 = ['K02635', 'K02637', 'K02634', 'K02636', 'K02642', 'K02643', 'K03689',
				'K02640']
#Cytochrome b6/f complex
	for i in cyt_b6:
		if i in ko_match:
			out_data['Cytochrome b6/f complex'] += 0.125
#Anoxygenic type-II reaction center pufL & pufM
	if ('K08928' in ko_match):
		out_data['anoxygenic type-II reaction center'] += 0.5
	if ('K08929' in ko_match):
		out_data['anoxygenic type-II reaction center'] += 0.5
#Anoxygenic type-I reaction center pscABCD
	rci = ['K08940', 'K08941', 'K08942', 'K08943']
	for i in rci:
		if i in ko_match:
			out_data['anoxygenic type-I reaction center'] += 0.25
#Retinal biosynthesis
	retinal = ['K06443', 'K02291', 'K10027', 'K13789']
	for i in retinal:
		if i in ko_match:
			out_data['Retinal biosynthesis'] += 0.25
#Retinal production from apo-carotenals
	if ('K00464' in ko_match):
		out_data['Retinal from apo-carotenals'] = 1
	return out_data




def entnerdoudoroff(ko_match):
	value = 0
#H6PD; hexose-6-phosphate dehydrogenase
	if 'K13937' in ko_match:
		total = 1
#edd; phosphogluconate dehydratase
		if 'K01690' in ko_match:
			total += 1
#2-dehydro-3-deoxyphosphogluconate aldolase
		if ('K01625' in ko_match or 'K17463' in ko_match or 'K11395' in ko_match):
			total += 1
		value = float(total) / float(3)
	else:
		total = 0
#G6PD; glucose-6-phosphate 1-dehydrogenase
		if 'K00036' in ko_match:
			total += 1
#6-phosphogluconolactonase
		if ('K01057' in ko_match or 'K07404' in ko_match):
			total += 1
#edd; phosphogluconate dehydratase
		if 'K01690' in ko_match:
			total += 1
#2-dehydro-3-deoxyphosphogluconate aldolase
		if ('K01625' in ko_match or 'K17463' in ko_match or 'K11395' in ko_match):
			total += 1
		value = float(total) / float(4)
	return {'Entner-Doudoroff Pathway': float("%.2f" % (value))}

def mixedacid(ko_match):
	out_data = {'Mixed acid: Lactate':0, 'Mixed acid: Formate':0,
	'Mixed acid: Formate to CO2 & H2':0, 'Mixed acid: Acetate':0,
	'Mixed acid: Ethanol, Acetate to Acetylaldehyde':0,
	'Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)':0,
	'Mixed acid: Ethanol, Acetylaldehyde to Ethanol':0,
	'Mixed acid: PEP to Succinate via OAA, malate & fumarate': 0}
#LDH; L-lactate dehydrogenase
	if 'K00016' in ko_match:
		out_data['Mixed acid: Lactate'] = 1
#pf1D; formate C-acetyltransferase
	if 'K00656' in ko_match:
		out_data['Mixed acid: Formate'] = 1
#formate dehydrogenase
	formatedh = ['K00122', 'K00125', 'K00126', 'K00123', 'K00124', 'K00127']
	for i in formatedh:
		if i in ko_match:
			out_data['Mixed acid: Formate to CO2 & H2'] += 0.167
#poxB; pyruvate dehydrogenase (quinone)
	if 'K00156' in ko_match:
		out_data['Mixed acid: Acetate'] = 1
#poxL; pyruvate oxidase + K01512 acyP; acylphosphatase
	if 'K00158' in ko_match:
		if out_data['Mixed acid: Acetate'] != 1:
			out_data['Mixed acid: Acetate'] += 0.5
#acyP; acylphosphatase
	if 'K01512' in ko_match:
		if out_data['Mixed acid: Acetate'] != 1:
			out_data['Mixed acid: Acetate'] += 0.5
#ACH1; acetyl-CoA hydrolase
	if 'K01067' in ko_match:
		out_data['Mixed acid: Acetate'] = 1
#pta; phosphate acetyltransferase or eutD; phosphotransacetylase
	if ('K13788' in ko_match or 'K04020' in ko_match):
		if (out_data['Mixed acid: Acetate'] != 1 and 'K00158' not in ko_match):
			out_data['Mixed acid: Acetate'] += 0.5
#lactate 2-monooxygenase
	if 'K00467' in ko_match:
		out_data['Mixed acid: Acetate'] = 1
#aldehyde dehydrogenase (NAD(P)+) OR aldB; aldehyde dehydrogenase
#aldehyde dehydrogenase (NAD+)
	aldehydedh = ['K00128', 'K14085', 'K00149', 'K00129', 'K00138']
	for i in aldehydedh:
		if i in ko_match:
			out_data['Mixed acid: Ethanol, Acetate to Acetylaldehyde'] = 1
#acetaldehyde dehydrogenase (acetylating)
	altaldehydedh = ['K00132', 'K04072', 'K04073', 'K18366', 'K04021']
	for i in altaldehydedh:
		if i in ko_match:
			out_data['Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)'] = 1
#alcohol dehydrogenase
	alchoholdh = ['K13951', 'K13980', 'K13952', 'K13953', 'K13954', 'K00001',
				'K00121', 'K04072', 'K18857', 'K00114', 'K00002', 'K04022']
	for i in alchoholdh:
		if i in ko_match:
			out_data['Mixed acid: Ethanol, Acetylaldehyde to Ethanol'] == 1
#methanol dehydrogenase
	if out_data['Mixed acid: Ethanol, Acetylaldehyde to Ethanol'] != 1:
		methanoldh = ['K14028', 'K14029']
		for i in methanoldh:
			if i in ko_match:
				out_data['Mixed acid: Ethanol, Acetylaldehyde to Ethanol'] += 0.5
#pckA; phosphoenolpyruvate carboxykinase (GTP) OR PEPCK; phosphoenolpyruvate carboxykinase (diphosphate) OR pckA; phosphoenolpyruvate carboxykinase (ATP)
	if ('K01596' in ko_match or 'K20370' in ko_match or 'K01610' in ko_match):
		out_data['Mixed acid: PEP to Succinate via OAA, malate & fumarate'] += 0.25
#malate dehydrogenase (NADP+) OR  mqo; malate dehydrogenase (quinone)
	if ('K00051' in ko_match or 'K00116' in ko_match):
		out_data['Mixed acid: PEP to Succinate via OAA, malate & fumarate'] += 0.25
#malate dehydrogenase
	if ('K00051' not in ko_match and 'K00116' not in ko_match):
		malatedh = ['K00025', 'K00026', 'K00024']
		for i in malatedh:
			if i in ko_match:
				out_data['Mixed acid: PEP to Succinate via OAA, malate & fumarate'] += 0.083
#fumarate hydratase, class I
	fumaratehydratase = ['K01676', 'K01677', 'K01678', 'K01679']
	for i in fumaratehydratase:
		if i in ko_match:
			out_data['Mixed acid: PEP to Succinate via OAA, malate & fumarate'] += 0.0625
#fumarate reductase flavoprotein
	fumaratereductase = ['K00244', 'K00245', 'K00246', 'K00247']
	for i in fumaratereductase:
		if i in ko_match:
			out_data['Mixed acid: PEP to Succinate via OAA, malate & fumarate'] += 0.0625
	return out_data

def naphthalene(ko_match):
	total = 0
#nahAabcd; naphthalene 1,2-dioxygenase
	nahAabcd = ['K14579', 'K14580', 'K14578', 'K14581']
	for i in nahAabcd:
		if i in ko_match:
			total += 0.25
#nahB; cis-1,2-dihydro-1,2-dihydroxynaphthalene/dibenzothiophene dihydrodiol dehydrogenase
	if 'K14582' in ko_match:
		total += 1
#nahC; 1,2-dihydroxynaphthalene dioxygenase
	if 'K14583' in ko_match:
		total += 1
#nahD; 2-hydroxychromene-2-carboxylate isomerase
	if 'K14584' in ko_match:
		total += 1
#nahE; trans-o-hydroxybenzylidenepyruvate hydratase-aldolase
	if 'K14585' in ko_match:
		total += 1
#nahF; salicylaldehyde dehydrogenase
	if 'K00152' in ko_match:
		total += 1
	value = float(total)/float(6)
	return {'Naphthalene degradation to salicylate': float("%.2f" % (value))}

def biofilm(ko_match):
	out_data = {'Biofilm PGA Synthesis protein':0,
	'Colanic acid and Biofilm transcriptional regulator':0,
	'Biofilm regulator BssS':0, 'Colanic acid and Biofilm protein A':0,
	'Curli fimbriae biosynthesis':0, 'Adhesion':0}
	pgasynth = ['K11935', 'K11931', 'K11936', 'K11937']
	for i in pgasynth:
		if i in ko_match:
			out_data['Biofilm PGA Synthesis protein'] += 0.25
	if 'K13654' in ko_match:
		out_data['Colanic acid and Biofilm transcriptional regulator'] = 1
	if 'K12148' in ko_match:
		out_data['Biofilm regulator BssS'] = 1
	if 'K13650' in ko_match:
		out_data['Colanic acid and Biofilm protein A'] = 1
	curli = ['K04335', 'K04334', 'K04336']
	for i in curli:
		if i in ko_match:
			out_data['Curli fimbriae biosynthesis'] += 0.33
	if 'K12687' in ko_match:
		out_data['Adhesion'] = 1
	return out_data

def competence(ko_match):
	out_data = {'Competence-related core components': 0,
				'Competence-related related components': 0,
				'Competence factors': 0}
	comp_core = ['K02237', 'K01493', 'K02238', 'K02239', 'K02240', 'K02241',
				'K02242', 'K02243', 'K02244', 'K02245', 'K02246', 'K02247',
				'K02248', 'K02249']
	for i in comp_core:
		if i in ko_match:
			out_data['Competence-related core components'] += 0.07
	comp_related = ['K02250', 'K02251', 'K02252', 'K02253', 'K02254']
	for i in comp_related:
		if i in ko_match:
			out_data['Competence-related related components'] += 0.2
	comp_factors = ['K12292', 'K07680', 'K12293', 'K12415', 'K12294',
					'K12295', 'K12296']
	for i in comp_factors:
		if i in ko_match:
			out_data['Competence factors'] += 0.14
	return out_data

def anaplerotic(ko_match):
	out_data = {'Glyoxylate shunt':0, 'Anaplerotic genes': 0}
#isocitrate lyase + malate synthase
	if 'K01637' in ko_match and 'K01638' in ko_match:
		out_data['Glyoxylate shunt'] = 1
#malate dehydrogenase (oxaloacetate-decarboxylating) (NADP+)
	if 'K00029' in ko_match:
		out_data['Anaplerotic genes'] += 0.25
#phosphoenolpyruvate carboxylase
	if 'K01595' in ko_match:
		out_data['Anaplerotic genes'] += 0.25
#phosphoenolpyruvate carboxykinase (ATP) or (GTP) or (diphosphate)
	if ('K01610' in ko_match) or ('K01596' in ko_match) or ('K20370' in ko_match):
		out_data['Anaplerotic genes'] += 0.25
#pyruvate carboxylase
	if ('K01958' in ko_match) or ('K01959' in ko_match and 'K01960' in ko_match):
		out_data['Anaplerotic genes'] += 0.25
	return out_data

def sulfolipid(ko_match):
	out_data = {'Sulfolipid biosynthesis':0}
	if 'K06118' in ko_match:
		out_data['Sulfolipid biosynthesis'] += 0.5
	if 'K06119' in ko_match:
		out_data['Sulfolipid biosynthesis'] += 0.5
	return out_data

def cplyase(ko_match):
	out_data = {'C-P lyase cleavage PhnJ':0, 'CP-lyase complex':0, 'CP-lyase operon':0}
#C-P lyase PhnJ
	if 'K06163' in ko_match:
		out_data['C-P lyase cleavage PhnJ'] = 1
#Tetradimer complex PhnJ, PhnG, PhnH, PhnI 
	complex_ = ['K06163', 'K06164', 'K06165', 'K06166']
	for i in complex_:
		if i in ko_match:
			out_data['CP-lyase complex'] += 0.25
#Full operon phnFGHIJKLMNOP - phosphonate transporter includes phnCED
	operon = ['K06163', 'K06164', 'K06165', 'K06166', 'K05780', 'K06162', 'K06167', 'K09994', 'K05774', 'K05781', 'K02043']
	for i in operon:
		if i in ko_match:
			out_data['CP-lyase operon'] += 0.09
	return out_data

def secretion(ko_match):
	out_data = {'Type I Secretion':0, 'Type III Secretion':0, 'Type II Secretion':0, 'Type IV Secretion':0, 'Type VI Secretion':0,
				'Sec-SRP':0, 'Twin Arginine Targeting':0, 'Type Vabc Secretion':0}
#Secretion mechanisms as described by KEGG Pathway Bacterial Secretion System
	typei = ['K12340', 'K11003', 'K11004']
	for i in typei:
		if i in ko_match:
			out_data['Type I Secretion'] += 0.33
	typeiii = ['K03221', 'K04056', 'K04057', 'K04059', 'K03219', 'K04058', 'K03222', 'K03226', 'K03227', 'K03228',
				'K03229', 'K03230', 'K03224', 'K03225', 'K03223']
	for i in typeiii:
		if i in ko_match:
			out_data['Type III Secretion'] += 0.0666
	typeii = ['K02453', 'K02465', 'K02452', 'K02455', 'K02456', 'K02457', 'K02458', 'K02459', 'K02460', 'K02461',
				'K02462', 'K02454', 'K02464']
	for i in typeii:
		if i in ko_match:
			out_data['Type II Secretion'] += 0.0769
	typeiv = ['K03194', 'K03197', 'K03198', 'K03200', 'K03202', 'K03204', 'K03201', 'K03203', 'K03195', 'K03199',
				'K03196', 'K03205']
	for i in typeiv:
		if i in ko_match:
			out_data['Type IV Secretion'] += 0.083
	typevi = ['K11904', 'K11903', 'K11906', 'K11891', 'K11892', 'K11907', 'K11912', 'K11913', 'K11915']
	for i in typevi:
		if i in ko_match:
			out_data['Type VI Secretion'] += 0.111
	tat = ['K03116', 'K03117', 'K03118', 'K03425']
	for i in tat:
		if i in ko_match:
			out_data['Twin Arginine Targeting'] += 0.25

	if ('K03072' in ko_match) or ('K03074' in ko_match):
		secsrp = ['K03072', 'K03074', 'K03073', 'K03075', 'K03076', 'K03210', 'K03217', 'K03070', 'K13301', 'K03110',
					'K03071', 'K03106']
		for i in secsrp:
			if i in ko_match:
				out_data['Sec-SRP'] += 0.083
	if ('K12257' in ko_match):
		secsrp = ['K12257', 'K03073', 'K03075', 'K03076', 'K03210', 'K03217', 'K03070', 'K13301', 'K03110',
					'K03071', 'K03106']
		out_data['Sec-SRP'] = 0
		for i in secsrp:
			if i in ko_match:
				out_data['Sec-SRP'] += 0.09
	typev = ['K11028', 'K11017', 'K11016', 'K12341', 'K12342']
	for i in typev:
		if i in ko_match:
			out_data['Type Vabc Secretion'] += 0.2
	return out_data

def serine(ko_match):
	out_data = {'Serine pathway/formaldehyde assimilation':0}
	serine_pathway = ['K00600', 'K00830', 'K00018', 'K11529', 'K01689', 'K01595',
						'K00024', 'K08692', 'K14067', 'K08692']
	for i in serine_pathway:
		if i in ko_match:
			out_data['Serine pathway/formaldehyde assimilation'] += .1
	return out_data

def arsenic(ko_match):
	out_data = {'Arsenic reduction':0}
#arsC
	if ('K00537' in ko_match) or ('K03741' in ko_match) or ('K18701' in ko_match):
		out_data['Arsenic reduction'] += 0.25
#arsB
	if ('K03325' in ko_match) or ('K03893' in ko_match):
		out_data['Arsenic reduction'] += 0.25
#arsR
	if 'K03892' in ko_match:
		out_data['Arsenic reduction'] += 0.25
#arsA
	if 'K01551' in ko_match:
		out_data['Arsenic reduction'] += 0.25
	return out_data

def metal_transport(ko_match):
	out_data = {'Cobalt transporter CbiMQ':0, 'Cobalt transporter CbtA':0,
	'Cobalt transporter CorA':0, 'Nickel ABC-type substrate-binding NikA':0,
	'Copper transporter CopA':0, 'Ferrous iron transporter FeoB': 0,
	'Ferric iron ABC-type substrate-binding AfuA': 0,
	'Fe-Mn transporter MntH':0}
#CbiMQ
	if 'K02007' in ko_match:
		out_data['Cobalt transporter CbiMQ'] += 0.5
	if 'K02008' in ko_match:
		out_data['Cobalt transporter CbiMQ'] += 0.5
#CbtA
	if 'K18837' in ko_match:
		out_data['Cobalt transporter CbtA'] = 1.0
#CorA
	if 'K03284' in ko_match:
		out_data['Cobalt transporter CorA'] = 1.0
#NikA
	if 'K15584' in ko_match:
		out_data['Nickel ABC-type substrate-binding NikA'] = 1.0
#CopA
	if 'K17686' in ko_match:
		out_data['Copper transporter CopA'] = 1.0
#FeoB
	if 'K04759' in ko_match:
		out_data['Ferrous iron transporter FeoB'] = 1.0
#AfuA
	if 'K02012' in ko_match:
		out_data['Ferric iron ABC-type substrate-binding AfuA'] = 1.0
#MntH
	if 'K03322' in ko_match:
		out_data['Fe-Mn transporter MntH'] = 1.0
	return out_data

def amino_acids(ko_match):
	out_data = {'histidine':0, 'arginine':0, 'lysine':0, 'serine':0, 
	'threonine':0, 'asparagine': 0, 'glutamine': 0, 'cysteine':0,
	'glycine':0, 'proline':0, 'alanine':0, 'valine':0, 
	'methionine':0, 'phenylalanine': 0, 'isoleucine': 0, 'leucine':0,
	'tryptophan':0, 'tyrosine': 0, 'aspartate': 0, 'glutamate':0}
	# histidine
	if 'K00013' in ko_match:
		out_data['histidine'] = 1
	# arginine
	if ('K01755' in ko_match) or ('K14681' in ko_match):
		out_data['arginine'] = 1
	# asparagine
	if ('K01913' in ko_match) or ('K01953' in ko_match):
		out_data['asparagine'] = 1
	# lysine
	lysine = ['K01586', 'K12526', 'K05831', 'K00290']
	for i in lysine:
		if i in ko_match:
			out_data['lysine'] = 1
	# serine
	serine = ['K01079', 'K02203', 'K02205']
	for i in serine:
		if i in ko_match:
			out_data['serine'] = 1
	if 'K00600' in ko_match:
		out_data['serine'] = 1
		out_data['glycine'] = 1
	# threonine
	if 'K01733' in ko_match:
		out_data['threonine'] = 1
	if 'K01620' in ko_match:
		out_data['threonine'] = 1
		out_data['glycine'] = 1
	# glutamine
	if 'K01915' in ko_match:
		out_data['glutamine'] = 1
	# cysteine
	cysteine = ['K01758', 'K17217', 'K01738', 'K10150', 'K12339']
	for i in cysteine:
		if i in ko_match:
			out_data['cysteine'] = 1
	# proline
	if ('K00286' in ko_match) or ('K01750' in ko_match):
		out_data['proline'] = 1
	# alanine
	if ('K14260' in ko_match) or ('K09758' in ko_match) or ('K00259' in ko_match) or ('K19244' in ko_match):
		out_data['alanine'] = 1
	# valine & isoleucine
	valine_isoleucine = ['K00826', 'K01687', 'K00053', 'K01652', 'K01653', 'K11258']
	for i in valine_isoleucine:
		if i in ko_match:
			out_data['valine'] += 0.166
			out_data['isoleucine'] += 0.166
	# leucine
	leucine = ['K00826', 'K00052', 'K01703', 'K01649']
	for i in leucine:
		if i in ko_match:
			out_data['leucine'] += 0.25
	# methionine
	if ('K00549' in ko_match) or ('K00548' in ko_match):
		out_data['methionine'] = 1
	# phenylalanine & tyrosine
	if ('K00832' in ko_match) or ('K00838' in ko_match):
		out_data['phenylalanine'] = 1
		out_data['tyrosine'] = 1
	#phenylalanine
	if ('K04518' in ko_match) or ('K05359' in ko_match) or ('K01713' in ko_match):
		out_data['phenylalanine'] = 1
	#tyrosine
	if ('K15226' in ko_match) or ('K24018' in ko_match) or ('K00220' in ko_match):
		out_data['tyrosine'] = 1
	# tryptophan
	if 'K01695' in ko_match:
		out_data['tryptophan'] += 0.5
	if ('K01696' in ko_match) or ('K06001' in ko_match):
		out_data['tryptophan'] += 0.5
	# aspartate & glutamate
	aspartate_glutamate = ['K00811', 'K00812', 'K00813', 'K11358', 'K14454', 'K14455']
	for i in aspartate_glutamate:
		if i in ko_match:
			out_data['aspartate'] = 1
			out_data['glutamate'] = 1

	return out_data

def plastic(ko_match):
	out_data = {"PET degradation": 0}
	#poly(ethylene terephthalate) hydrolase
	#mono(ethylene terephthalate) hydrolase
	#1,2-dihydroxy-3,5-cyclohexadiene-1,4-dicarboxylate dehydrogenase
	petdeg = ["K21104", "K21105", "K18076"]
	for i in petdeg:
		if i in ko_match:
			out_data["PET degradation"] += 0.25
	#terephthalate 1,2-dioxygenase oxygenase
	#two possible versions
	if ("K18077" in ko_match) or ("K18074" in ko_match and "K18075" in ko_match):
		out_data["PET degradation"] += 0.25

	return out_data

def carbon_storage(ko_match):
    out_data = {'starch/glycogen synthesis': 0, 'starch/glycogen degradation': 0, 'polyhydroxybutyrate synthesis': 0}
    # starch synthesis
    carbonsto = ["K00703", "K00975"]
    for i in carbonsto:
        if i in ko_match:
            out_data["starch/glycogen synthesis"] += 0.33

    if ('K00700' in ko_match) or ('K16149' in ko_match):
        out_data['starch/glycogen synthesis'] += 0.33
    #starch > D-glucose
		#K21574	susB; glucan 1,4-alpha-glucosidase
    if ('K21574' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    #starch > cyclodextrin
    #K00701	cgt; cyclomaltodextrin glucanotransferase
    if ('K00701' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    #starch > maltodextrin
    #K01214	treX; isoamylase
    if ('K01214' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    #starch > glucose-6P
    if ('K00688' in ko_match) or ('K16153' in ko_match) or ('K00705' in ko_match) or ('K22451' in ko_match) or ('K02438' in ko_match) or ('K01200' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    #starch > dextrin
    #alpha-amlyase
    if ('K01176' in ko_match) or ('K05343' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    #beta-amlyase
    if 'K01177' in ko_match:
        out_data['starch/glycogen degradation'] = 1
    #maltogenic alpha-amylase
    if ('K05992' in ko_match) or ('K01208' in ko_match):
        out_data['starch/glycogen degradation'] = 1
    if ('K00023' in ko_match):
        out_data['polyhydroxybutyrate synthesis'] += 0.5
    phb = ['K00626', 'K03821', 'K22881']
    for i in phb:
        if i in ko_match:
            out_data['polyhydroxybutyrate synthesis'] += 0.167

    return out_data


def phosphate_storage(ko_match):
    out_data = {'bidirectional polyphosphate': 0}

    if ('K00937' in ko_match) or ('K22468' in ko_match):
        out_data['bidirectional polyphosphate'] += 0.5
    if ('K01507' in ko_match) or ('K15986' in ko_match) or ('K06019' in ko_match):
        out_data['bidirectional polyphosphate'] += 0.5

    return out_data

def carotenoids(ko_match):
	out_data = {'carotenoids backbone biosynthesis':0, 'end-product astaxanthin':0,
	'end-product nostoxanthin':0, 'end-product zeaxanthin diglucoside':0, 
	'end-product myxoxanthophylls': 0, 'staphyloaxanthin biosynthesis':0,
	'mevalonate pathway': 0, 'MEP-DOXP pathway':0}

	staphyloaxanthin = ['K10208', 'K10209', 'K10210', 'K00128',
	'K10211', 'K10212']
	for i in staphyloaxanthin:
		if i in ko_match:
			out_data['staphyloaxanthin biosynthesis'] += 0.167
	mevalonate = ['K01641', 'K00054', 'K00869','K00938', 'K01597']
	for i in mevalonate:
		if i in ko_match:
			out_data['mevalonate pathway'] += 0.2
	mepdoxp = ['K01662', 'K00099', 'K00991', 'K00919', 'K01770',
	'K03526', 'K03527']
	for i in mepdoxp:
		if i in ko_match:
			out_data['MEP-DOXP pathway'] += 0.142
	backbone = ['K23155', 'K02291', 'K02291', 'K00514']
	for i in backbone:
		if i in ko_match:
			out_data['carotenoids backbone biosynthesis'] += 0.2
	if ('K06443' in ko_match) or ('K14606' in ko_match):
		out_data['carotenoids backbone biosynthesis'] += 0.2
	if ('K02294' in ko_match):
		out_data['end-product nostoxanthin'] = 1
	if ('K02294' in ko_match) and ('K14596' in ko_match):
		out_data['end-product zeaxanthin diglucoside'] = 1
	if ('K09836' in ko_match) or ('K02292' in ko_match):
		out_data['end-product astaxanthin'] = 1
	myxoxanthophylls = ['K08977', 'K02294', 'K00721']
	for i in myxoxanthophylls:
		if i in ko_match:
			out_data['end-product myxoxanthophylls'] += 0.33

	return out_data



def default_viz(genome_df, outfile_name):
	import seaborn as sns
	import matplotlib.pyplot as plt
	sns.set(font_scale=1.2)
	sns.set_style({"savefig.dpi": 200})
	ax = sns.heatmap(genome_df, cmap=plt.cm.YlOrRd, linewidths=2,
		linecolor='k', square=True, xticklabels=True,
		yticklabels=True, cbar_kws={"shrink": 0.1})
	ax.xaxis.tick_top()
	#ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
	plt.xticks(rotation=90)
	plt.yticks(rotation=0)
	# get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
	fig = ax.get_figure()
	# specify dimensions and save
	#xLen = len(genome_df.columns.values.tolist())*20
	#yLen = len(genome_df.index.tolist())*20
	fig.set_size_inches(100, 100)
	fig.savefig(outfile_name, bbox_inches='tight', pad_inches=0.1)

def main():
	import os
	import matplotlib
	matplotlib.use('Agg')
	import argparse
	import pandas as pd
	from scipy.cluster import hierarchy
	from scipy.spatial import distance


	parser = argparse.ArgumentParser(description="Accepts KEGG KOALA\
									text file as input. Produces function\
									list and heat map figure.")
	parser.add_argument('-i', '--input', help="Input KOALA file. See documentation\
						for correct format")
	parser.add_argument('-t', '--tangleopt', help="Number of tree iterations for minimizing tangles in tanglegram", default=1000)
	parser.add_argument('-o', '--output', help="List version of the final heat\
						map figure")
	parser.add_argument('-v', '--vizoption', help="Options: static, interactive, tanglegram")
	parser.add_argument('--newick', help="Required input for tanglegram visualization")
	parser.add_argument("-m", "--myorder", help ="Orders output as specified by	user.", default="None")
	args = parser.parse_args()
	arg_dict = vars(args)

	genome_data = {}

	for line in open(str(arg_dict['input']), "r"):
		line = line.rstrip()
		info = line.split()
		if len(info) > 1:
			if info[0].split("_")[0] in genome_data.keys():
				genome_data[info[0].split("_")[0]].append(info[1])
			else:
				genome_data[info[0].split("_")[0]] = [info[1]]

	function_order = ['glycolysis', 'gluconeogenesis', 'TCA Cycle',
	'NAD(P)H-quinone oxidoreductase', 'NADH-quinone oxidoreductase',
	'Na-NADH-ubiquinone oxidoreductase',
	'F-type ATPase', 'V-type ATPase', 'Cytochrome c oxidase',
	'Ubiquinol-cytochrome c reductase', 'Cytochrome o ubiquinol oxidase',
	'Cytochrome aa3-600 menaquinol oxidase',
	'Cytochrome c oxidase, cbb3-type', 'Cytochrome bd complex', 'RuBisCo',
	'CBB Cycle', 'rTCA Cycle', 'Wood-Ljungdahl',
	'3-Hydroxypropionate Bicycle', '4-Hydroxybutyrate/3-hydroxypropionate',
	'pectinesterase', 'diacetylchitobiose deacetylase', 'glucoamylase',
	'D-galacturonate epimerase', 'exo-poly-alpha-galacturonosidase',
	'oligogalacturonide lyase', 'cellulase', 'exopolygalacturonase',
	'chitinase', 'basic endochitinase B', 'bifunctional chitinase/lysozyme',
	'beta-N-acetylhexosaminidase', 'D-galacturonate isomerase',
	'alpha-amylase', 'beta-glucosidase', 'pullulanase',
	'ammonia oxidation (amo/pmmo)', 'hydroxylamine oxidation', 'nitrite oxidation',
	'dissim nitrate reduction', 'DNRA', 'nitrite reduction',
	'nitric oxide reduction', 'nitrous-oxide reduction',
	'nitrogen fixation', 'hydrazine dehydrogenase', 'hydrazine synthase',
	'dissimilatory sulfate < > APS',
	'dissimilatory sulfite < > APS', 'dissimilatory sulfite < > sulfide',
	'thiosulfate oxidation', 'alt thiosulfate oxidation tsdA',
	'alt thiosulfate oxidation doxAD', 'sulfur reductase sreABC',
	'thiosulfate/polysulfide reductase', 'sulfhydrogenase',
	'sulfur disproportionation', 'sulfur dioxygenase',
	'sulfite dehydrogenase', 'sulfite dehydrogenase (quinone)',
	'sulfide oxidation', 'sulfur assimilation',
	'DMSP demethylation', 'DMS dehydrogenase', 'DMSO reductase',
	'NiFe hydrogenase', 'ferredoxin hydrogenase',
	'membrane-bound hydrogenase', 'hydrogen:quinone oxidoreductase', 'NAD-reducing hydrogenase',
	'NADP-reducing hydrogenase', 'NiFe hydrogenase Hyd-1',
	'thiamin biosynthesis',
	'riboflavin biosynthesis' ,
	'cobalamin biosynthesis', 'transporter: vitamin B12',
	'transporter: thiamin', 'transporter: urea',
	'transporter: phosphonate', 'transporter: phosphate',
	'Flagellum', 'Chemotaxis', 'Methanogenesis via methanol',
	'Methanogenesis via acetate',
	'Methanogenesis via dimethylsulfide, methanethiol, methylpropanoate',
	'Methanogenesis via methylamine', 'Methanogenesis via trimethylamine',
	'Methanogenesis via dimethylamine', 'Methanogenesis via CO2',
	'Coenzyme B/Coenzyme M regeneration',
	'Coenzyme M reduction to methane', 'Soluble methane monooxygenase',
	'methanol dehydrogenase', 'alcohol oxidase',
	'dimethylamine/trimethylamine dehydrogenase',
	'Photosystem II', 'Photosystem I', 'Cytochrome b6/f complex',
	'anoxygenic type-II reaction center', 'anoxygenic type-I reaction center',
	'Retinal biosynthesis', 'Retinal from apo-carotenals',
	'Entner-Doudoroff Pathway', 'Mixed acid: Lactate', 'Mixed acid: Formate',
	'Mixed acid: Formate to CO2 & H2', 'Mixed acid: Acetate',
	'Mixed acid: Ethanol, Acetate to Acetylaldehyde',
	'Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)',
	'Mixed acid: Ethanol, Acetylaldehyde to Ethanol',
	'Mixed acid: PEP to Succinate via OAA, malate & fumarate',
	'Naphthalene degradation to salicylate',
	'Biofilm PGA Synthesis protein',
	'Colanic acid and Biofilm transcriptional regulator',
	'Biofilm regulator BssS', 'Colanic acid and Biofilm protein A',
	'Curli fimbriae biosynthesis', 'Adhesion', 'Competence-related core components',
	'Competence-related related components', 'Competence factors',
	'Glyoxylate shunt', 'Anaplerotic genes', 'Sulfolipid biosynthesis',
	'C-P lyase cleavage PhnJ', 'CP-lyase complex', 'CP-lyase operon', 'Type I Secretion',
	'Type III Secretion', 'Type II Secretion', 'Type IV Secretion', 'Type VI Secretion',
	'Sec-SRP', 'Twin Arginine Targeting', 'Type Vabc Secretion',
	'Serine pathway/formaldehyde assimilation', 'Arsenic reduction',
	'Cobalt transporter CbiMQ', 'Cobalt transporter CbtA',
	'Cobalt transporter CorA', 'Nickel ABC-type substrate-binding NikA',
	'Copper transporter CopA', 'Ferrous iron transporter FeoB',
	'Ferric iron ABC-type substrate-binding AfuA',
	'Fe-Mn transporter MntH', 'histidine', 'arginine', 'lysine', 'serine', 
	'threonine', 'asparagine', 'glutamine', 'cysteine',
	'glycine', 'proline', 'alanine', 'valine', 
	'methionine', 'phenylalanine', 'isoleucine', 'leucine',
	'tryptophan', 'tyrosine', 'aspartate', 'glutamate', 'PET degradation',
	'starch/glycogen synthesis', 'starch/glycogen degradation', 'polyhydroxybutyrate synthesis',
	'bidirectional polyphosphate', 'carotenoids backbone biosynthesis', 
	'end-product astaxanthin', 'end-product nostoxanthin', 'end-product zeaxanthin diglucoside', 
	'end-product myxoxanthophylls', 'staphyloaxanthin biosynthesis',
	'mevalonate pathway', 'MEP-DOXP pathway']


	filehandle = str(arg_dict['output'])
	out_file = open(filehandle, "w")
	out_file.write('Function'+"\t"+str("\t".join(function_order))+"\n")

	for k in genome_data:
		pathway_data = {}
		pathway_data.update(nitrogen(genome_data[k]))
		pathway_data.update(glycolysis(genome_data[k]))
		pathway_data.update(gluconeogenesis(genome_data[k]))
		pathway_data.update(tca_cycle(genome_data[k]))
		pathway_data.update(cbb_cycle(genome_data[k]))
		pathway_data.update(reverse_tca(genome_data[k]))
		pathway_data.update(wood_ljungdahl(genome_data[k]))
		pathway_data.update(three_prop(genome_data[k]))
		pathway_data.update(four_hydrox(genome_data[k]))
		pathway_data.update(c_degradation(genome_data[k]))
		pathway_data.update(chemotaxis(genome_data[k]))
		pathway_data.update(flagellum(genome_data[k]))
		pathway_data.update(sulfur(genome_data[k]))
		pathway_data.update(methanogenesis(genome_data[k]))
		pathway_data.update(methane_ox(genome_data[k]))
		pathway_data.update(hydrogen(genome_data[k]))
		pathway_data.update(transporters(genome_data[k]))
		pathway_data.update(riboflavin(genome_data[k]))
		pathway_data.update(thiamin(genome_data[k]))
		pathway_data.update(oxidative_phoshorylation(genome_data[k]))
	#Addendum 2
		pathway_data.update(photosynthesis(genome_data[k]))
		pathway_data.update(entnerdoudoroff(genome_data[k]))
		pathway_data.update(mixedacid(genome_data[k]))
		pathway_data.update(naphthalene(genome_data[k]))
		pathway_data.update(biofilm(genome_data[k]))
		pathway_data.update(cobalamin(genome_data[k]))
		pathway_data.update(competence(genome_data[k]))
		pathway_data.update(anaplerotic(genome_data[k]))
		pathway_data.update(sulfolipid(genome_data[k]))
		pathway_data.update(cplyase(genome_data[k]))
		pathway_data.update(secretion(genome_data[k]))
		pathway_data.update(serine(genome_data[k]))
		pathway_data.update(arsenic(genome_data[k]))
		pathway_data.update(metal_transport(genome_data[k]))
		pathway_data.update(amino_acids(genome_data[k]))
		pathway_data.update(plastic(genome_data[k]))
		pathway_data.update(carbon_storage(genome_data[k]))
		pathway_data.update(phosphate_storage(genome_data[k]))
		pathway_data.update(carotenoids(genome_data[k]))
	#    print k, pathway_data

		out_string = str(k)+"\t"
		out_list = [k]
		for i in function_order:
			out_list.append(pathway_data[i])
		out_string = str(out_list).strip('[]')
		tab_string = ""
		for l in out_string:
			if l == "\'":
				continue
			if l == ",":
				tab_string = tab_string + "\t"
			else:
				tab_string = tab_string + l
		out_file.write(tab_string+"\n")
	out_file.close()


	file_in = open(filehandle, "r")
	genome = pd.read_csv(file_in, index_col=0, sep='\t')
	rearrange = False
	if arg_dict["myorder"] != 'None' and os.path.exists(arg_dict["myorder"]):
		rearrange = True
		leaf_order = []
		for line in open(str(arg_dict["myorder"]), "r"):
			line = line.rstrip("\r\n")
			leaf_order.append(line)
		genome = genome.reindex(leaf_order)

	if arg_dict['vizoption'] == 'static':
		from .KEGG_clustering import hClust_euclidean
		#from KEGG_clustering import hClust_euclidean
		if len(genome.index) >= 2 and not rearrange:
			genome = hClust_euclidean(genome)
		default_viz(genome, os.path.splitext(filehandle)[0] + ".svg")
	if arg_dict['vizoption'] == 'interactive':
		from .Plotly_viz import plotly_viz
		#from Plotly_viz import plotly_viz
		plotly_viz(genome, os.path.splitext(filehandle)[0] + ".html")
	if arg_dict['vizoption'] == 'tanglegram':
		from .MakeTanglegram import make_tanglegram
		#from MakeTanglegram import make_tanglegram
		if len(genome.index) >= 3:
			make_tanglegram(genome, str(arg_dict['newick']), os.path.splitext(filehandle)[0] + ".tanglegram.svg", int(arg_dict["tangleopt"]))
		else:
			raise ValueError("Tanglegram mode requires three or more genomes")


if __name__ == "__main__":
	main()
