#name: Chemical Space Using UMAP
#description: Chemical space using Uniform Manifold Approximation and Projection
#help-url: https://datagrok.ai/help/domains/chem/functions/umap
#language: python
#sample: chem/smiles_coordinates.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
#input: int neighbors = 15 [Number of neighbours]
#input: int minClusterSize = 3 [Minimum cluster size]
#output: graphics spanningTree
#output: graphics linkageTree
#output: graphics chemSpace

import umap
import numba
import hdbscan
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs

smiles = data[smiles]
mols = [Chem.MolFromSmiles(mol) for mol in smiles]
mols = [mol for mol in mols if mol is not None]
for mol in mols:
    AllChem.Compute2DCoords(mol)

X = []
for mol in mols:
    arr = np.zeros((0,))
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    DataStructs.ConvertToNumpyArray(fp, arr)
    X.append(arr)

@numba.njit()
def tanimoto_dist(a, b):
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc

umap_X = umap.UMAP(n_neighbors=neighbors, min_dist=0.3, metric=tanimoto_dist).fit_transform(X)
cluster_umap = hdbscan.HDBSCAN(min_cluster_size=minClusterSize, gen_min_span_tree=True)
cluster_umap.fit(umap_X)

plt.figure(0)
cluster_umap.minimum_spanning_tree_.plot(
    edge_cmap='viridis', edge_alpha=0.6,
    node_size=90, edge_linewidth=2)

plt.figure(1)
cluster_umap.single_linkage_tree_.plot(cmap='viridis', colorbar=True)

plt.figure(2)
plt.scatter(umap_X.T[0], umap_X.T[1], c=cluster_umap.labels_, cmap='plasma')
plt.title('Chemical space')
