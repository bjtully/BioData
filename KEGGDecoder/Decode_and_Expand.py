#!/usr/bin/python

'''
Combines the *.list format output files of KEGG-decoder.py and
KEGG-expander.py to form a single figure of all functions


'''

import argparse

parser = argparse.ArgumentParser(description="Accepts HMM search results of expander_dbvX.hmm\
								text file as input. Produces function\
								list and heat map figure.")
parser.add_argument('KOALA LIST', help="Input KOALA function list format. As generated from\
					KEGG-decoder")
parser.add_argument('HMM LIST', help="Input HMM function list format. As generated from\
					KEGG-expander")
args = parser.parse_args()
arg_dict = vars(args)

import matplotlib.pyplot as plt

import pandas as pd


koala = pd.read_table(open(str(arg_dict['KOALA LIST']), "r"), index_col=0)
hmm = pd.read_table(open(str(arg_dict['HMM LIST']), "r"), index_col=0)
output_df = koala.merge(hmm, left_index=True, right_index=True)

#Reorganize column orientation to put like pathways together
cols = output_df.columns.tolist()
retinal_index = cols.index('Retinal biosynthesis')
cols.insert(retinal_index+1, cols.pop(int(cols.index('beta-carotene 15,15-monooxygenase'))))
cols.insert(retinal_index+2, cols.pop(int(cols.index('rhodopsin'))))
trans_urea = cols.index('transporter: urea')
cols.insert(trans_urea+1, cols.pop(int(cols.index('transporter: ammonia'))))
nifH_index = cols.index('nitrogen fixation')
cols.insert(nifH_index+1, cols.pop(int(cols.index('Vanadium-only nitrogenase'))))
cols.insert(nifH_index+2, cols.pop(int(cols.index('Iron-only nitrogenase'))))
dmsplyase_index = cols.index('DMSP demethylation')
cols.insert(dmsplyase_index, cols.pop(int(cols.index('DMSP lyase (dddLQPDKW)'))))
cols.insert(dmsplyase_index+1, cols.pop(int(cols.index('DMSP synthase (dsyB)'))))
output_df = output_df[cols]

import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style({"savefig.dpi": 200})
ax = sns.heatmap(output_df, cmap=plt.cm.YlOrRd, linewidths=2, linecolor='k', square=True)
ax.xaxis.tick_top()
#ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
# get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
fig = ax.get_figure()
# specify dimensions and save
fig.set_size_inches(100, 100)
fig.savefig("decode-expand_heatmap.svg")