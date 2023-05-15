---
title: Chemical functions
---

## Supported functions

| Category         | Name                               | Function                  |
|------------------|------------------------------------|---------------------------|
|Cheminformatics   | Substructure search                |<br /><pre>`\#{x.ChemSubstructureSearch}`</pre>|
|Cheminformatics   | [PLACEHOLDER]                      |<br /><pre>`\#{x.ChemFindMCS}`</pre> |
|Cheminformatics   | Descriptors                        |<br /><pre>`\#{x.ChemDescriptors}`</pre> |
|Cheminformatics   | R-Groups                           |<br /><pre>`\#{x.ChemGetRGroups}`</pre> |
|Cheminformatics   | Fingerprints                       |<br /><pre>`\#{x.ChemFingerprints}`</pre> |
|Cheminformatics   | [PLACEHOLDER]                      |<br /><pre>`\#{x.ChemSimilaritySPE}`</pre> |
|Cheminformatics   | SMILES to InchI                    |<br /><pre>`\#{x.ChemSmilesToInchi}`</pre> |
|Cheminformatics   | SMILES to Canonical                |<br /><pre>`\#{x.ChemSmilesToCanonical}`</pre> |
|Cheminformatics   | Chemical map identifiers           |<br /><pre>`\#{x.ChemMapIdentifiers}`</pre> |
|Chemical analysis | Butina cluster                     |<br /><pre>`\#{x.ChemScripts:ButinaMoleculesClustering}`</pre> |
|Chemical analysis | Filter by catalogs                 |<br /><pre>`\#{x.ChemScripts:FilterByCatalogs}`</pre> |
|Chemical analysis | Gasteiger partial charges          |<br /><pre>`\#{x.ChemScripts:GasteigerPartialCharges}`</pre> |
|Chemical analysis | Murcko scaffolds                   |<br /><pre>`\#{x.ChemScripts:MurckoScaffolds}`</pre>|
|Chemical analysis | Similarity maps using fingerprints |<br /><pre>`\#{x.ChemScripts:SimilarityMapsUsingFingerprints}`</pre> |
|Chemical analysis | Chemical space using tSNE          |<br /><pre>`\#{x.ChemScripts:ChemicalSpaceUsingtSNE}`</pre> |
|Chemical analysis | Reactions                          |<br /><pre>`\#{x.ChemScripts:TwoComponentReaction}`</pre> |
|Chemical analysis | Chemical space using UMAP          |<br /><pre>`\#{x.ChemScripts:ChemicalSpaceUsingUMAP}`</pre> |
|Chemical analysis | USRCAT                             |<br /><pre>`\#{x.ChemScripts:USRCAT}`</pre> |
|TBD               | Mutate                             |<br /><pre>`[PLACEHOLDER]`</pre> |
|TBD               | Solubility prediction              |<br /><pre>`\#{x.18b704d0-0b50-11e9-b846-1fa94a4da5d1."Predict Solubility"}`</pre>|

### Butina cluster

Uses desired similarity within the cluster, as defined by Tanimoto index, as the only input to the clustering program.

References:

* [RDKit Cluster Butina Module](https://rdkit.org/docs/source/rdkit.ML.Cluster.Butina.html)
* [Butina JCICS 39 747-750 (1999)](http://www.l4patterns.com/uploads/dbclus-paper.pdf)

### Chemical space using tSNE

tSNE, short for t-distributed Stochastic Neighbor Embedding, is a data visualization tool designed to handle high-dimensional data. It achieves this by transforming the similarities between data points into joint probabilities, then minimizing the Kullback-Leibler divergence between the low-dimensional embedding and the original high-dimensional data. tSNE uses a non-convex cost function, meaning that different initializations can lead to different results. The following image illustrates the use of tSNE to visualize chemical space.

![Chemical Space Using tSNE](../../../uploads/chem/tsne.png "Chemical Space Using tSNE")

References:

* [RDKit](https://www.rdkit.org)
* [tSNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)

### Chemical space using UMAP

Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for
visualisation similarly to [tSNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html), but
also for general non-linear dimension reduction.

![Chemical Space Using UMAP](../../../uploads/chem/umap.png "Chemical Space Using UMAP")

References:

* [RDKit](https://www.rdkit.org)
* [UMAP](https://umap-learn.readthedocs.io/en/latest/)
* [tSNE scikit](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)

### Filter by catalogs

Screen out or reject undesireable molecules based on various criteria.

Filter sets:

* **PAINS**: Pan assay interference patterns, seperated into three sets (PAINS_A, PAINS_B, and PAINS_C).
* **BRENK**: Filters unwanted functionality due to potential toxicity reasons or unfavorable pharmacokinetics.
* **NIH**: Annotated compounds with problematic functional groups
* **ZINC**: Filtering based on drug-likeness and unwanted functional groups.

References:

* [RDKit FilterCatalogs](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FilterCatalog/README)

### Gasteiger partial charges

Visualizes atomic charges in a molecule.

![Gasteiger Partial Charges](../../../uploads/chem/gasteiger-charges.png "Gasteiger Partial Charges")

References:

* [RDKit Visualization of Descriptors](https://www.rdkit.org/docs/GettingStartedInPython.html#visualization-of-descriptors)
* [Gasteiger-Marsili empirical atomic partial charges](https://www.codessa-pro.com/descriptors/electrostatic/gmc.htm)

### Murcko scaffolds

Converts a column with molecules to Murcko scaffolds.

![Murcko Scaffolds](../../../uploads/chem/murcko-scaffolds.png "Murcko Scaffolds")

References:

* [rdkit.Chem.Scaffolds.MurckoScaffold module](https://rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html)
* [Computational Exploration of Molecular Scaffolds in Medicinal Chemistry](https://europepmc.org/abstract/MED/26840095)
* [Comparative analyses of structural features and scaffold diversity for purchasable compound libraries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5400773/)

### Mutate

Mutate molecules using different mechanisms:

* Adding atoms
* Adding bonds
* Removing bonds

Mutations can be randomized using **randomize** flag. Mutation mechanisms and place will be in randomized for each mutation step.

References:

* [RDKit](https://www.rdkit.org/)

### Reactions

Reaction template is in [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) format. Reactants can be combined from two sets, or sequentially depending on the **matrixExpansion** flag.

![Reactions](../../../uploads/chem/reactions.png "Reactions")

References:

* [RDKit Chemical reaction handling](https://rdkit.org/docs/RDKit_Book.html#chemical-reaction-handling)
* [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)

### Similarity maps using fingerptints

Visualizes the atomic contributions to the similarity between a molecule and a reference molecule.

![Similarity Maps Using Fingerprints](../../../uploads/chem/sim-maps.png "Similarity Maps Using Fingerprints")

References:

* [RDKit generating similarity maps using fingerprints](https://www.rdkit.org/docs/GettingStartedInPython.html#generating-similarity-maps-using-fingerprints)
* [Similarity maps - a visualization strategy for molecular fingerprints and machine-learning methods](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43)

### Solubility prediction

The H2O modeling engine was used to train the model using the "Solubility Train" dataset <br />(`#{x.Demo:SolubilityTrain."Solubility Train"}`). The modelling method used was "Generalized Linear Modeling".

Molecular descriptors used in the model:

* **MolWt**: Molecular weight
* **Ipc**: The information content of the coefficients of the characteristic polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule
* **TPSA**: Total polar surface area
* **LabuteASA**: Labute's approximate surface area
* **NumHDonors**: Number of hydrogen donors
* **NumHAcceptors**: Number of hydrogen acceptors
* **MolLogP**: Wildman-Crippen LogP value
* **HeavyAtomCount**: Number of heavy atoms
* **NumRotatableBonds**: Number of rotatable bonds
* **RingCount**: Number of rings
* **NumValenceElectrons**: Number of valence electrons

References:

* [RDKit](https://www.rdkit.org)
* [MoleculeNet: A benchmark for molecular machine learning](https://arxiv.org/abs/1703.00564)

### USRCAT

USRCAT is an extension of the Ultrafast Shape Recognition (USR) algorithm, which is used for molecular shape-based virtual screening to discover new chemical scaffolds in compound libraries. USRCAT incorporates pharmacophoric information in addition to molecular shape, which enables it to distinguish between compounds with similar shapes but distinct pharmacophoric features.

![USRCAT](../../../uploads/chem/usrcat.png "USRCAT")

References:

* [RDKit](https://www.rdkit.org)
* [USRCAT](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3505738/)
