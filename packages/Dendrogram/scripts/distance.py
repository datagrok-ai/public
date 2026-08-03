#name: distanceScript
#friendlyName: Distance Matrix
#description: Compute a condensed pairwise distance matrix between rows of the input table.
#language: python
#input: dataframe data [Input data table]
#input: string distance_name = 'euclidean' {choices: ['euclidean', 'manhattan']} [Distance metric]
#output: dataframe result

import numpy as np
import pandas as pd

from scipy.spatial.distance import cdist, pdist
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree

column_array = data[data.columns].to_numpy()
dist_list = pdist(column_array, distance_name)

result = pd.DataFrame.from_dict({'distance': dist_list})
