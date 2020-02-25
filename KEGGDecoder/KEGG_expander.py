#!/usr/bin/env python3

'''
KEGG-expander.py V.0.5
V.0.5 Removal of Amphibactin biosynthesis componenets. Additions of dsrD to help distinguish
oxidative and reductive dsrAB
V.0.4.1 Adds argument that should allow script to generate figure on most Unix systems
V.0.4 Adds amphibactin biosynthesis, ferrioxamine biosynthesis
V.0.3.1 Added parameters to force labels to be printed.
KEGG-expander.py V.0.3
V.0.3. Adds DMSP lyase (dddQ, dddP, dddD, dddK, dddW), DMSP synthase (dsyB)
Usage: python KEGG-decoder.py <HMM TBL INPUT> <FUNCTION LIST FORMAT>

Designed to parse through the hmmsearch results table generated from
the expander_dbvX.hmm to generate a heatmap figure similar to
KEGG-decoder.py 

Recommended to run hmmsearch as follows:
hmmsearch --tblout <HMM TBL NAME> -T 75 expander_dbvX.hmm <PROTEIN SEQS>
bit score cutoff of 75, equivalent to e-value < 10^-20

Dependencies:
Pandas - http://pandas.pydata.org/pandas-docs/stable/install.html
Seaborn - http://seaborn.pydata.org/installing.html
matplotlib - http://matplotlib.org/users/installing.html

For extended information about HMM assignments, genes and pathways,
please see accompanying document "Pfam_definitions.txt"


'''


def main():
	import matplotlib
	matplotlib.use('Agg')
	import argparse

	parser = argparse.ArgumentParser(description="Accepts HMM search results of expander_dbvX.hmm\
									text file as input. Produces function\
									list and heat map figure.")
	parser.add_argument('Input', help="Input HMM table file. See documentation\
						for correct format")
	parser.add_argument('Output', help="List version of the final heat\
						map figure")
	args = parser.parse_args()
	arg_dict = vars(args)

	genome_data = {}

	for line in open(str(arg_dict['Input']), "r"):
		if line[0] != "#":
			line = line.rstrip()
			info = line.split()
			genome_id = info[0].split("_")[0]
			# Sfams for DMSP lyase dddP and dddD require a more strigent bit score cutoff (>500)
			if info[3].split(".")[0] == "14591" or info[3].split(".")[0] == "25993":
				if float(info[5]) > 500:
					try:
						genome_data[genome_id].append(info[3].split(".")[0])
					except KeyError:
						genome_data[genome_id] = [info[3].split(".")[0]]
				else:
					continue
			if info[3].split(".")[0] == "4254":
				if float(info[5]) > 260:
					try:
						genome_data[genome_id].append(info[3].split(".")[0])
					except KeyError:
						genome_data[genome_id] = [info[3].split(".")[0]]
				else:
					continue
			else:
				try:
					genome_data[genome_id].append(info[3].split(".")[0])
				except KeyError:
					genome_data[genome_id] = [info[3].split(".")[0]]
			# Sfams for amphibactin biosynthesis requires a more stringent bit score cutoff (>1000)
			#		if info[3].split(".")[0] == "1544" or info[3].split(".")[0] == "27549":
			#			if float(info[5]) > 1000:
			#				try:
			#					genome_data[genome_id].append(info[3].split(".")[0])
			#				except KeyError:
			#					genome_data[genome_id] = [info[3].split(".")[0]]
			#			else:
			#				continue
			# Sfams for ferrioxamine biosynthesis requires a more stringent bit score cutoff (>200)
			if info[3].split(".")[0] == "2219" or info[3].split(".")[0] == "2732" or info[3].split(".")[0] == "9429" or \
					info[3].split(".")[0] == "51934":
				if float(info[5]) > 200:
					try:
						genome_data[genome_id].append(info[3].split(".")[0])
					except KeyError:
						genome_data[genome_id] = [info[3].split(".")[0]]
				else:
					continue

	def rhodopsin(hmm_match):
		out_data = {'beta-carotene 15,15-monooxygenase': 0, 'rhodopsin': 0}
		if 'PF01036' in hmm_match:
			out_data['rhodopsin'] = 1
		if 'TIGR03753' in hmm_match:
			out_data['beta-carotene 15,15-monooxygenase'] = 1
		return out_data

	# def peptidases(hmm_match):
	#	out_data = {'Peptidase family C25': 0, 'Bacterial pre-peptidase C-terminal domain': 0,
	#		'Clostripain family': 0, 'Peptidase family M28': 0, 'Peptidase family M50': 0,
	#		'Di- and tripeptidases': 0, 'Leucyl aminopeptidase': 0, 'Xaa-Pro aminopeptidase': 0,
	#		'Peptidase propeptide and YPEB domain': 0, 'Oligoendopeptidase F': 0,
	#		'Phosphoserine aminotransferase': 0, 'Lipoprotein signal peptidase': 0,
	#		'Aminopeptidase N': 0, 'Zinc carboxypeptidase': 0, 'Peptidase S24-like': 0,
	#		'Peptidase S26': 0, 'D-aminopeptidase': 0, 'M61 glycyl aminopeptidase': 0}
	#	if 'PF01364' in hmm_match:
	#		out_data['Peptidase family C25'] = 1
	#	if 'PF04151' in hmm_match:
	#		out_data['Bacterial pre-peptidase C-terminal domain'] = 1
	#	if 'PF03415' in hmm_match:
	#		out_data['Clostripain family'] = 1
	#	if 'PF04389' in hmm_match:
	#		out_data['Peptidase family M28'] = 1
	#	if 'PF02163' in hmm_match:
	#		out_data['Peptidase family M50'] = 1
	#	if 'PF01546' in hmm_match:
	#		out_data['Di- and tripeptidases'] = 1
	#	if 'PF02073' in hmm_match:
	#		out_data['Leucyl aminopeptidase'] = 1
	#	if 'PF00557' in hmm_match:
	#		out_data['Xaa-Pro aminopeptidase'] = 1
	#	if 'PF03413' in hmm_match:
	#		out_data['Peptidase propeptide and YPEB domain'] = 1
	#	if 'PF01432' in hmm_match:
	#		out_data['Oligoendopeptidase F'] = 1
	#	if 'PF00266' in hmm_match:
	#		out_data['Phosphoserine aminotransferase'] = 1
	#	if 'PF01252' in hmm_match:
	#		out_data['Lipoprotein signal peptidase'] = 1
	#	if 'PF01433' in hmm_match:
	#		out_data['Aminopeptidase N'] = 1
	#	if 'PF00246' in hmm_match:
	#		out_data['Zinc carboxypeptidase'] = 1
	#	if 'PF00717' in hmm_match:
	#		out_data['Peptidase S24-like'] = 1
	#	if 'PF10502' in hmm_match:
	#		out_data['Peptidase S26'] = 1
	#	if 'PF04951' in hmm_match:
	#		out_data['D-aminopeptidase'] = 1
	#	if 'PF05299' in hmm_match:
	#		out_data['M61 glycyl aminopeptidase'] = 1
	#	return out_data

	def alt_nitrogenase(hmm_match):
		out_data = {'Vanadium-only nitrogenase': 0, 'Iron-only nitrogenase': 0}
		v_nitro = ['TIGR01860', 'TIGR02932', 'TIGR02930']
		for i in v_nitro:
			if i in hmm_match:
				out_data['Vanadium-only nitrogenase'] += 0.33
		fe_nitro = ['TIGR01861', 'TIGR02931', 'TIGR02929']
		for i in fe_nitro:
			if i in hmm_match:
				out_data['Iron-only nitrogenase'] += 0.33
		return out_data

	def amm_trans(hmm_match):
		out_data = {'transporter: ammonia': 0}
		if 'PF00909' in hmm_match:
			out_data['transporter: ammonia'] = 1
		return out_data

	def dmsplyase(hmm_match):
		out_data = {'DMSP lyase (dddLQPDKW)': 0}
		dmsp = ['PF16867', '14591', '25993', '94923', '274874']
		for i in dmsp:
			if i in hmm_match:
				out_data['DMSP lyase (dddLQPDKW)'] = 1
		return out_data

	def dmspsynthase(hmm_match):
		out_data = {'DMSP synthase (dsyB)': 0}
		if '4254' in hmm_match:
			out_data['DMSP synthase (dsyB)'] = 1
		return out_data

	# def amphibactin(hmm_match):
	#	out_data = {'amphibactin ACO2092-3homolog':0}
	#	if ('1544' in hmm_match) and ('27549' in hmm_match):
	#		out_data['amphibactin ACO2092-3homolog'] = 1
	#	return out_data

	def ferrioxamine(hmm_match):
		out_data = {'ferrioxamine biosynthesis': 0}
		ferrioxamine = ["2219", "2732", "9429", "51934"]
		for i in ferrioxamine:
			if i in hmm_match:
				out_data['ferrioxamine biosynthesis'] += 0.25
		return out_data

	def dissim_sulfite(hmm_match):
		out_data = {'DsrD dissimilatory sulfite reductase': 0}
		if 'PF08679' in hmm_match:
			out_data['DsrD dissimilatory sulfite reductase'] = 1
		return out_data

	def metal_transport(hmm_match):
		out_data = {'Cobalt transporter CbtB': 0, 'Copper binding HMA protein': 0,
					'Fe-Zn-Mn permease ZupT': 0}
		if 'TIGR02459' in hmm_match:
			out_data['Cobalt transporter CbtB'] = 1
		if 'TIGR00003' in hmm_match:
			out_data['Copper binding HMA protein'] = 1
		if 'TIGR00820' in hmm_match:
			out_data['Fe-Zn-Mn permease ZupT'] = 1

		return out_data

	function_order = ['beta-carotene 15,15-monooxygenase', 'rhodopsin', 'Vanadium-only nitrogenase',
					  'Iron-only nitrogenase', 'transporter: ammonia',
					  'DMSP lyase (dddLQPDKW)', 'DMSP synthase (dsyB)',
					  'ferrioxamine biosynthesis', 'DsrD dissimilatory sulfite reductase', 'Cobalt transporter CbtB',
					  'Copper binding HMA protein', 'Fe-Zn-Mn permease ZupT']

	filehandle = str(arg_dict['Output'])
	out_file = open(filehandle, "w")
	out_file.write('Function' + "\t" + str("\t".join(function_order)) + "\n")

	for k in genome_data:
		pathway_data = {}
		pathway_data.update(rhodopsin(genome_data[k]))
		#	pathway_data.update(peptidases(genome_data[k]))
		pathway_data.update(alt_nitrogenase(genome_data[k]))
		pathway_data.update(amm_trans(genome_data[k]))
		pathway_data.update(dmsplyase(genome_data[k]))
		pathway_data.update(dmspsynthase(genome_data[k]))
		#	pathway_data.update(amphibactin(genome_data[k]))
		pathway_data.update(ferrioxamine(genome_data[k]))
		pathway_data.update(dissim_sulfite(genome_data[k]))
		pathway_data.update(metal_transport(genome_data[k]))

		out_string = str(k) + "\t"
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
		out_file.write(tab_string + "\n")

	out_file.close()

	import matplotlib.pyplot as plt

	import pandas as pd

	file_in = open(filehandle, "r")
	genome = pd.read_csv(file_in, index_col=0, sep='\t')
	import seaborn as sns
	sns.set(font_scale=1.2)
	sns.set_style({"savefig.dpi": 200})
	ax = sns.heatmap(genome, cmap=plt.cm.YlOrRd, linewidths=2, linecolor='k', square=True, xticklabels=True,
					 yticklabels=True)
	ax.xaxis.tick_top()
	# ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
	plt.xticks(rotation=90)
	plt.yticks(rotation=0)
	# get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
	fig = ax.get_figure()
	# specify dimensions and save
	fig.set_size_inches(100, 100)
	fig.savefig("hmm_heatmap.svg")


if __name__ == "__main__":
	main()
