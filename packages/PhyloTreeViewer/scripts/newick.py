#name: getNewick
#description: Returns the newick representation of the tree for given dataset
#language: python
#input: dataframe data [Input data table]
#output: string newick

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

column_array = data[data.columns].to_numpy()
dist_matrix = cdist(column_array, column_array)
one_dimension = pdist(dist_matrix, 'euclidean')
link_matrix = linkage(one_dimension, 'single')
leaf_names = list(range(0, len(column_array)))
tree = to_tree(link_matrix, False)
newick = get_newick(tree, tree.dist, leaf_names)