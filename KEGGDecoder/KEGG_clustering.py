#!/usr/bin/python

'''
This script is a Heirarchical clustering module for the KEGG-decoder.py
versions after V.0.8
Runs scipy clustering with various metrix on the KEGG_Decoder genome DataFrame
'''
from scipy.cluster.hierarchy import ward, complete, average, dendrogram, fcluster, linkage

def hClust_euclidean(genome_df):
	linkage_matrix = linkage(genome_df, method='average', metric='euclidean')
	#linkage_matrix = linkage(df, metric='braycurtis')
	names = genome_df.index.tolist()
	#clust = dendrogram(linkage_matrix, orientation="right", labels=names, get_leaves=True)
	clust = dendrogram(linkage_matrix, no_plot=True, labels=names, get_leaves=True)
	leaves = clust['ivl']
	leave_order = list(leaves)
	genome_df = genome_df.reindex(leave_order)

	return genome_df

def hClust_correlation(genome_df):
	linkage_matrix = linkage(genome_df, method='single', metric='correlation')
	#linkage_matrix = linkage(df, metric='braycurtis')
	names = genome_df.index.tolist()
	#clust = dendrogram(linkage_matrix, orientation="right", labels=names, get_leaves=True)
	clust = dendrogram(linkage_matrix, no_plot=True, labels=names, get_leaves=True)
	leaves = clust['ivl']
	leave_order = list(leaves)
	genome_df = genome_df.reindex(leave_order)

	return genome_df

def hClust_most_least(genome_df):
	sort_dex = genome_df.sum(axis=1).sort_values(ascending=True).index
	genome_df = genome_df.ix[sort_dex]

	return genome_df

def hClust_least_most(genome_df):
	sort_dex = genome_df.sum(axis=1).sort_values(ascending=False).index
	genome_df = genome_df.ix[sort_dex]

	return genome_df

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description="This file is intended as a Plotly module for the KEGG_decoder")
	args = parser.parse_args()
	arg_dict = vars(args)