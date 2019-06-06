#!/usr/bin/python

'''
This script is an implementation for tanglegram for the KEGG-decoder.py
versions after V.0.8
Runs tanglegram on two DataFrames - one generated from clustered KEGG
metabolisms & one phylogenetic newick file provided by the user
Added by Taylor Reiter : tereiter@ucdavis.edu
'''

def make_tanglegram(genome_df, newick):
    import matplotlib.pyplot as plt
    import itertools
    from Bio import Phylo
    import tanglegram as tg
    from scipy.spatial.distance import pdist, squareform
    
    # FORMAT KEGGDECODER OUTPUT
    # generate distance matrix for genome_df from pathway values
    genome_df = pd.read_csv(genome_df, index_col=0, sep='\t')
    kegg_d = squareform(pdist(genome_df, metric='euclidean'))
    kegg_m = pd.DataFrame(kegg_d)
    kegg_m.columns = genome_df.index.tolist()
    kegg_m.index = genome_df.index.tolist()
    kegg_m = kegg_m.reindex(sorted(kegg_m.columns), axis=1) # reorder column names alphabetically
    kegg_m.sort_index(inplace=True) # reorder row names alphabetically
    
    # FORMAT NEWICK FILE
    # generate distance matrix from newick file
    tree = Phylo.read(newick, 'newick')
    
    tree_d = {}
    for x, y in itertools.combinations(tree.get_terminals(), 2):
        v = tree.distance(x, y)
        tree_d[x.name] = tree_d.get(x.name, {})
        tree_d[x.name][y.name] = v
        tree_d[y.name] = tree_d.get(y.name, {})
        tree_d[y.name][x.name] = v
        
    for x in tree.get_terminals():
        tree_d[x.name][x.name] = 0
    
    tree_m = pd.DataFrame(tree_d)
    tree_m = tree_m.reindex(sorted(tree_m.columns), axis=1) # reorder column names alphabetically
    tree_m.sort_index(inplace=True) # reorder row names alphabetically
    
    # TANGLEGRAM
    kegg_labels = kegg_m.columns.values.tolist()
    tree_labels = tree_m.columns.values.tolist()

    kegg_mat = pd.DataFrame(kegg_m,
                        columns=kegg_labels,
                        index=kegg_labels)

    tree_mat = pd.DataFrame(tree_m,
                        columns=tree_labels,
                        index=tree_labels)

    # Plot and try to minimize cross-over
    fig = tg.gen_tangle(kegg_mat, tree_mat, optimize_order=1000)
    fig.set_size_inches(10, 10)
    fig.savefig("function_newick_tanglegram.svg")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="This file is intended as a Tanglegram module for the KEGG_decoder")
    args = parser.parse_args()
    arg_dict = vars(args)