<!-- TITLE: Cheminformatics -->

# Cheminformatics

Datagrok provides a powerful set of tools for [Cheminformatics](https://en.wikipedia.org/wiki/Cheminformatics), accelerating chemically-related workflows.

## Overview

With Datagrok, you can:

* Instantly import your data of the most common [molecular structure formats](../../access/supported-formats.md#molecular-structure-formats).
  For databases and cloud services we provide [30+ data connectors](../../access/data-connection.md), [data queries](../../access/data-query.md) and [data preparation pipelines](../../access/data-pipeline.md).  
  Our [semantic types](../../discover/semantic-types.md) include molecules, enabling automatic SMILES to structures rendering, chemically-aware viewers and adaptive info panels.

  #GIF с демонстрацией формата SMILES + рендеринг молекул + изменяющиеся инфо панели

* Visualize your data with context-driven [Viewers](../../visualize/viewers/). For molecular data, [chemically aware viewers](chemically-aware-viewers) offer additional functionality:
  * Vizualization of the molecules on grid, scatter plot, histogram, trellis plot

    #GIF, показывающая, как выглядят структуры на осях графиков>

  * Synchronized filtering across selected viewers.
  
    #GIF с одновременной фильтрацией на всех вьюерах>
  
  * Datapoint associated tooltip with useful information
  
    #GIF с наведением на датапоинты>

* [Sketcher](sketcher.md) lets you sketch a new molecule, edit an existing one, or retrieve one by entering compound identifier/trivial name.

  ![Sketcher](../../uploads/chem/sketcher.png "Sketcher")

  While drawing a molecule, info panels update interactively.
  #GIF с обновляющимися инфо-панелями>

  With a structure in Sketcher, substructure search is available. You can select or filter rows in dataset by sketched substuctures.

  #GIF с процессом поиска субструктур в молекуле>
  
  We support several Sketcher types: Marvin, ChemDraw, OpenChemLib, Ketcher.

* Analyze molecules in [info panels](../../discover/info-panels.md): view structures in 2D and 3D, evaluate [drug likeness](info-panels/drug-likeness.md), [toxicity risks](info-panels/toxicity-risks), [structural alerts](info-panels/structural-alerts.md), chemical properties

* Explore molecules with [diversity search](diversity-search.md) and [similarity search](similarity-search.md)
* Extract [molecular descriptors](descriptors.md) and [fingerprints](fingerprints.md) on the fly
* Train models using [predictive modeling](chem-predictive-modeling.md) and incorporate them in pipelines or [info panels](../../discover/info-panels.md)

<!-- ![Cell](cell-renderer_test.gif "Automatic SMILES rendering") -->

## Augment

Observe new information from your data in [info panels](../../discover/info-panels.md). Datagrok's info panels collect molecule-specific functions in one place, calculating them on the fly.


In info panels all the functions fall into Chemistry, Biology, Databases groups.

<tab>

* Identifiers - fetches all known identifiers for the specified structure across [UniChem databases](https://www.ebi.ac.uk/unichem/)
* Structure 2D – visualizes a molecule in 2D
* Structure 3D – visualizes a molecule in 3D through generating .mol file
* Molfile – provides a .mol file
* [Gasteiger Partial Charges](functions/gasteiger-charges.md) - visualizes and highlights partial charges in a molecule
* [Chem descriptors](descriptors.md) – calculates and displays the specified descriptors for a molecule
* Properties – yields the list of calculated or predicted physical and chemical properties: empirical formula, molecular weight, hydrogen bond acceptor (HBA) and donor (HDA) values, LogP and LogS, polar surface area (PSA), number of rotatable bonds and stereocenters, also IUPAC name.
* [Toxicity](info-panels/toxicity-risks.md) –  drug design related feature to predict the toxicity scores. Consists of such categories as mutagenicity, tumorogenicity, irritating effects, reproductive effects.
* [Structural alerts](info-panels/structural-alerts.md) – drug design related feature to highlight the fragments in structures that might greatly increase the toxicity and other problematic structural features
* [Drug likeness](info-panels/structural-alerts.md) – drug design related feature to get a score that shows how likely this molecule is to be a drug.

  ![Toxicity, Gasteiger Partial Charges, Solubility Prediction](../../uploads/gifs/chem-model-augment.gif "Toxicity, Gasteiger Partial Charges, Solubility Prediction")

In addition to these predefined info panels, users also can develop their own using any scripting language supported by the Datagrok platform.

## Calculate

With chemical dataset uploaded, you can conduct necessary calculations to move along your workflow. In Datagrok, two types of calculations are available: descriptors/fingerprints and mapping functions.

### Descriptors and fingerprints

To vectorize molecular graph data, the platform supports [descriptors](descriptors.md) and [fingerprints](fingerprints.md).

<!-- With fingerprints and descriptors calculated, the following tasks become available: -->
While descriptors are more physical and fingerprints are more abstract vectors, both enable the following procedures:

* Similarity and diversity search
* Chem space dimensionality reduction
* SAR analysis
* Machine learning predictive modeling

Among supported descriptors are Lipinski, Crippen, EState, EState VSA, Fragments, Graph, MolSurf, QED and [others](descriptors.md).

#GIF, высчитывающий fingerprints !!! подождать, пока фингерпринты зальют в top-menu Chem>

Among supported fingerprints are RDKFingerprint, MACCSKeys, AtomPair, TopologicalTorsion, Morgan/Circular and [others](fingerprints.md).

![descriptors](descriptors.gif "Calculate descriptors")

### Mapping functions

Effective management across chemistry data is available only when the same structures are stored under the same identifiers (when the same structures are interlinked across databases)

For each molecule, mapping functions calculate its unique textual identifiers known as [International Chemical Identifiers](https://en.wikipedia.org/wiki/International_Chemical_Identifier).

With just a few clicks, you can convert your structures to InChI and InChI keys, its hashed version.

#GIF c конверсиями>

## Transform

Datagrok lets transform your molecular structures for standartization and/or augmentation purposes.
Beside standard [data wrangling](../../transform/data-wrangling.md) procedures, you can use chemistry-specific ones: curation and mutation.

### Curation

Variant 1:
[Сhemical structure curation](chem-curate.md) standardizes your chemical structures and thereby improves your SAR analysis or your model's prediction accuracy.
Moreover, by chemical structure curation you can avoid data-associated errors such as duplicated vectors in the training set or incorrect structure representation.

Variant 2:
[Сhemical structure curation](chem-curate.md) standardizes your chemical structures and thereby helps avoid data-associated errors such as duplicated vectors in the training set or incorrect structure representation.
By chemical structure curation you can improve your SAR analysis or your model's prediction accuracy.

We offer the following curation methods:

* Kekulization
* Reionization
* Normalization
* Neutralization
* Tautomerization
* Main fragment selection

![Chemical structure curation](chem_curation_demo.gif)

### Mutation

Mutation function generates new structures based on the specified one. Generation is combinatorial and includes such options as step - number of permutations applied to the structure, and max random result - number of output structures.

<!-- пока про Mutation не стоит писать, ибо функция сыровата -->
<br>
## Search

Оur search functions help to navigate across molecular structures and include the following types:

* ### [Similarity and diversity search](similarity-search.md)

  For the specified structure, finds the 10 most similar/diverse ones in your dataset.
  In similarity search all the found structures are sorted by similarity score  based on Morgan fingerprints and multiple [distance metrics](similarity-spe#available-distance-metrics) to opt.

  ![Similarity search](../../uploads/gifs/similarity-search.gif "Similarity search")

* ### [Substructure search](substructure-search.md)

  Sketcher-supported function that searches a specified structural pattern in the datasource.

  In Datagrok, substructure search implements advanced navigation with unique features:

  * Sketcher-oriented structural adjustment of the found molecules
  * Explicit hydrogen search support
  * Extended toolset for aromaticity search
  * Multicolumn filtering

> Similarity and substructure search support uploaded dataset, public databases (Chembl, PubChem, DrugBank) and [relational databases](db-substructure-similarity-search.md).

#GIF с поиском субструктур в опен-сорсных БД

> Note: to expand search capabilities, use Sketcher to find sketched substructures in an uploaded dataset or in public databases.

#GIF c поиском с помощью скетчера>

## Analyze

Var 1:
The core of Grok's cheminformatics domain is our analytical toolset. For chemistry-related tasks we offer multiple methods of 2 functional categories:

Var 2:
To handle chemical analysis, Grok's toolset offers multiple methods. Functionally, they fall into 2 categories:

Var 3:
Analytical methods vary. Thus, for easy navigation our analytical toolset falls into two functional categories:

* Analyze SAR - identify bioactivity values based on structure
* Analyze Structure - investigate exclusively structural pecularities

<!-- Leonid: 
* Analyze SAR - methods identifing structure-activity relationshpis
* Analyze Structure - approaches investigating structural pecularities of your dataset -->

### Analyze SAR

* #### [R-group analysis](r-group-analysis.md)

  Identifies R-group branches connected to a scaffold. For comparative visualizations, trellis plot is available.

  ![gt1](../../uploads/chem/graphTools2.png "rGroups") - сделать гифку

* Activity Cliffs

  Finds similar compounds with different activity in your dataset.

  Choose similarity percent of neighbors, run the function and explore activity cliffs displayed on 2D plots.

  #GIF с процедурой поиска Activity cliffs

* Structural Alerts

  Flags potential chemical hazards taken from several rule sets containing
[1,251 substructures in the SMARTS format](https://raw.githubusercontent.com/PatWalters/rd_filters/5f70235b387baa39669f25d95079e5dfec49a47c/rd_filters/data/alert_collection.csv).
  Presence of any of these substructures triggers a structural alert shown on the info panel.

### Analyze Structure

Analyze Structure methods include:

* Chem space

  Explores structural similarity of your molecular data on 2D clusters. For dimensionality reduction, you can choose either t-SNE or UMAP algorithm. Among distance metrics are Cosine, Sokal, Tanimoto, Assimetric.

* Elemental Analysis
  
  Represents molecules as the sum of its constituent atoms with optional radar plots.

  > Use elemental analysis while filtering molecules with inappropriate atoms

* Scaffold Tree

  Constructs tree hierarchy based on molecular scaffolds in your dataset. Use the hierarchy to filter or select the corresponding rows in the dataset.
  
## Predict

Datagrok   enables   machine   learning   predictive   models   by   using
chemical properties, descriptors, and fingerprints as features, and the
observed properties as targets. Within Datagrok, you can train, assess, execute,
share and use your models in pipelines.

We offer several models: PCA, Naive Bayes, K-Means, Distributed Random Forest, Gradient Boosting, XGBoost, Generalized Linear Modeling, AutoML and Deep Learning.

#Показать, как вставляются нововычисленные значения в info panels