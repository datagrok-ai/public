<!-- TITLE: Tests: Jupyter Notebook -->
<!-- SUBTITLE: -->

# Tests: Jupyter Notebook

Jupyter Notebook is an open-source web application that allows you to create and share documents that contain live code,
equations, visualizations and narrative text.

## Jupyter Notebooks

1. Open "smiles_coordinates.csv" table

1. Select "smiles_coordinates.csv" table on [Property Panel](../overview/navigation.md#properties)

1. Click "Open in Jupyter" in "Actions" tab on [Property Panel](../overview/navigation.md#properties)

* "smiles_coordinates.csv" table open in Jupyter
* Jupyter Notebook open in new tab

1. Click on "Open in new window" button on Jupyter Notebook view

* Current notebook open new browser window

1. Return to platform browser window

1. Click on "Open as HTML" button on Jupyter Notebook view

* Current notebook open in new tab in HTML markup

1. Return to Jupyter Notebook view

1. Enter script to "In" field in Jupyter Notebook (*)

1. Click "Run" in Jupyter Notebook

* After script execution, result is displayed in "Out\[1\]" field
* Result contains three graphical outputs

1. Click "Save" button in Jupyter Notebook

* Current notebook saved

1. Open [Workspace](../overview/workspace.md) window

* Notebook entity is in the "Scratchpad" project

1. Upload project with "smiles_coordinates.csv" table and created Notebook

1. Close everything open in platform

1. Go to "Datasets..." page on "Welcome" view and open and open project from step 12

* Table view and Notebook view are open

1. Go to Notebook view tab

* Notebook contains script saved in step 10 with result "Out\[1\] field

1. Edit current Notebook and save

1. Repeat steps 12-15 to check Notebook save changes

## Notebooks browser

1. Open table "smiles_cordinates"

1. Open **Tools | Notebooks | Browse Notebooks**

* New view "Notebooks" added
* "Card" is default view for Notebooks Browser
* There are "Name", "Description", "Author", "Created" and "Tags" on notebook's cards
* If suitable table is open to use in notebook,this is displayed on it's card

1. Change view to "Brief"

* Notebooks are displayed briefly

1. Call context menu for "Chemical Space Using tSNE" notebook

* Actions "Open", "Edit", "Delete", "Save as JSON", "Rename" are available from context menu
* Notebook can be applied to "smiles_cordinates" table from context menu

1. Choose "Chemical Space Using tSNE" notebook for [Properties Panel](../overview/navigation.md#properties)

* Tabs "Details", "Projects", "Activity", "Actions" and "Chats" are present
  on [Properties Panel](../overview/navigation.md#properties)

1. Click on "Created by me" filter on [Toolbox](../overview/navigation.md#toolbox)

* Only notebooks created by current user are displayed in browser
* Search field filled with predefined query *author = @current*

1. Test other filters

1. Change sort order

(*):

```import datagrok.api as dapi
   grok_host = 'https://dev.datagrok.ai/api'
   data = dapi.readTable(grok_host, '3ce28310-3260-11e9-b75e-71dd1cd6b072')
   %matplotlib inline
   %config InlineBackend.figure_format = 'retina'
   import matplotlib
   matplotlib.rcParams['figure.figsize'] = (8.0, 8.0)
   smiles = "smiles"
   components = 2
   minClusterSize = 3
   import numba
   import hdbscan
   import numpy as np
   import matplotlib.pyplot as plt
   from sklearn.manifold import TSNE
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

   tsne = TSNE(n_components=components, metric=tanimoto_dist)
   tsne_X = tsne.fit_transform(X)
   cluster_tsne = hdbscan.HDBSCAN(min_cluster_size=minClusterSize, gen_min_span_tree=True)
   cluster_tsne.fit(tsne_X)

   plt.figure(0)
   cluster_tsne.minimum_spanning_tree_.plot(
       edge_cmap='viridis', edge_alpha=0.6,
       node_size=90, edge_linewidth=2)

   plt.figure(1)
   cluster_tsne.single_linkage_tree_.plot(cmap='viridis', colorbar=True)

   plt.figure(2)
   plt.scatter(tsne_X.T[0], tsne_X.T[1], c=cluster_tsne.labels_, cmap='plasma')
   plt.title('Chemical space')```
