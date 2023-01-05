#name: hierarchicalClusteringByDistanceScript
#description: Returns the newick representation of the tree for given distance matrix
#language: python
#input: dataframe data [Input distance matrix condensed]
#input: double size [Input size (obs number)]
#input: string linkage_name = 'ward' {choices: ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']} [Linkage]
#output: string newick

import sys
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

one_dimension = np.array(data.loc[:, 'distance'])
link_matrix = linkage(one_dimension, linkage_name)
leaf_names = list(range(0, size))
tree = to_tree(link_matrix, False)
newick = get_newick(tree, tree.dist, leaf_names)