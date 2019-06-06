#!/usr/bin/python

'''
This script is a Plotly vizualization module written for KEGG-decoder.py 
versions after V.0.8 Builds an interactive heatmap of metabolic pathways
specified by KEGG-decoder.py
Added by Roth Conrad : rotheconrad@gatech.edu
'''

def plotly_viz(genome_df):
	# build heatmap in plotly.offline
	from .KEGG_clustering import hClust_euclidean
	Euclidean_genome_df = hClust_euclidean(genome_df)

	from .KEGG_clustering import hClust_correlation
	Correlation_genome_df = hClust_correlation(genome_df)

	from .KEGG_clustering import hClust_most_least
	Most_Least_genome_df = hClust_most_least(genome_df)

	from .KEGG_clustering import hClust_least_most
	Least_Most_genome_df = hClust_least_most(genome_df)

	import plotly.graph_objs as go
	import plotly.offline as py

	xLen = len(genome_df.columns.values.tolist())*20
	yLen = len(genome_df.index.tolist())*20

	colorscale = [
					[0, '#f1eef6'],
					[0.2, '#f1eef6'],
					[0.2 ,'#bdc9e1'],
					[0.4 ,'#bdc9e1'],
					[0.4 ,'#74a9cf'],
					[0.6 ,'#74a9cf'],
					[0.6 ,'#2b8cbe'],
					[0.8 ,'#2b8cbe'],
					[0.8 ,'#045a8d'],
					[1 ,'#045a8d']]

	colorbar = {'tick0': 0, 'dtick': 0.2, 'lenmode': 'pixels', 'len': 500, 'y': 1}

	Euclidean_clust = go.Heatmap(x=Euclidean_genome_df.columns.values.tolist(), 
						y=Euclidean_genome_df.index.tolist(), 
						z=Euclidean_genome_df.values.tolist(), 
						colorscale=colorscale, 
						colorbar=colorbar, 
						xgap = 1, 
						ygap = 1)

	Correlation_clust = go.Heatmap(x=Correlation_genome_df.columns.values.tolist(), 
						y=Correlation_genome_df.index.tolist(), 
						z=Correlation_genome_df.values.tolist(), 
						colorscale=colorscale, 
						colorbar=colorbar, 
						xgap = 1, 
						ygap = 1,
						visible=False)

	Most_Least_clust = go.Heatmap(x=Most_Least_genome_df.columns.values.tolist(), 
						y=Most_Least_genome_df.index.tolist(), 
						z=Most_Least_genome_df.values.tolist(), 
						colorscale=colorscale, 
						colorbar=colorbar, 
						xgap = 1, 
						ygap = 1,
						visible=False)

	Least_Most_clust = go.Heatmap(x=Least_Most_genome_df.columns.values.tolist(), 
						y=Least_Most_genome_df.index.tolist(), 
						z=Least_Most_genome_df.values.tolist(), 
						colorscale=colorscale, 
						colorbar=colorbar, 
						xgap = 1, 
						ygap = 1,
						visible=False)

	data = [Euclidean_clust, Correlation_clust, Most_Least_clust, Least_Most_clust]

	updatemenus = [dict(
					buttons = [
					dict(label = 'Euclidean_Clustering', method = 'update', args = [{'visible': [True, False, False, False]}]),
					dict(label =  'Correlation_Clustering', method = 'update', args = [{'visible': [False, True, False, False]}]),
					dict(label =  'Most_to_Least', method = 'update', args = [{'visible': [False, False, True, False]}]),
					dict(label =  'Least_to_Most', method = 'update', args = [{'visible': [False, False, False, True]}])
					], 
					direction = 'down',
					pad = {'r': 10, 't': 10},
					showactive = True,
					x = 0.1,
					xanchor = 'left',
					y = 1.1,
					yanchor = 'top'
					)]

	layout = go.Layout(xaxis={'side': 'top'},
						autosize=False,
						width=xLen,
						height=yLen,
						plot_bgcolor='#000000',
						margin=go.layout.Margin(t=500),
						updatemenus=updatemenus,
						)



	fig = go.Figure(data=data, layout=layout)
	py.plot(fig, filename='function_heatmap.html')
	# py.iplot(data, filename='pandas.heatmap')

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="This file is intended as a Plotly module for the KEGG_decoder")
	args = parser.parse_args()
	arg_dict = vars(args)
