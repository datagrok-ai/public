---
title: "Cheminformatics"
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Datagrok provides a powerful set of tools for [Cheminformatics](https://en.wikipedia.org/wiki/Cheminformatics), accelerating chemically-related workflows.

## Overview

With Datagrok, you can:

* Instantly import your data in [multiple data formats](../../access/supported-formats.md), including the most common [molecular structure formats](../../access/supported-formats.md#molecular-structure-formats).  
  To access databases and cloud services, you can use [30+ data connectors](../../access/data-connection.md), [data queries](../../access/data-query.md), and [data preparation pipelines](../../access/data-pipeline.md).  

  By recognizing molecules in your data, the platform automatically renders chemical structures and connects adaptive [info panels](../../discover/info-panels.md).

  ![SMILES rendering](smiles-rendering.gif "SMILES rendering")

* Visualize your data with [viewers](../../visualize/viewers/). Our viewers, including [chemically aware viewers](chemically-aware-viewers.md) for molecular data, offer the following functionality:

  * Visualization of molecules on axes.

    ![Molecules on axis](structures-on-axis.gif "Molecules on axis")

  * Synchronized filtering across selected viewers.

    ![Synchronized filtering](synchronized-filtering.gif "Synchronized filtering")

  * Useful tooltips for each datapoint.

    ![Datapoint tooltip](datapoint-tooltip.gif "Datapoint tooltip")

* Sketch a new molecule, edit an existing one, or retrieve one by entering compound identifier/trivial name with [sketcher](sketcher.md). Datagrok supports the following sketchers: Marvin, ChemDraw, OpenChemLib, Ketcher.

  ![Sketcher](sketcher-molecules.gif "Molecular structures in sketcher")
  
  <details>
  <summary>Use sketcher to filter your data by a substructure.</summary>
  <div>

  ![Sketcher substructure search](sketcher-substructure-search.gif "Sketcher substructure search")

  </div>
  </details>

* Explore molecules in [info panels](../../discover/info-panels.md): view structures in 2D and 3D, evaluate [drug likeness](info-panels/drug-likeness.md), [toxicity risks](info-panels/toxicity-risks), [structural alerts](info-panels/structural-alerts.md), chemical properties.

* Search through a collection of molecules with [diversity search](diversity-search.md), [similarity search](similarity-search.md), and substructure search.

* Calculate [molecular descriptors](descriptors.md) and [fingerprints](fingerprints.md).

* Train models using [predictive modeling](chem-predictive-modeling.md) and incorporate them in pipelines or [info panels](../../discover/info-panels.md).

## Augment

Extract more information about your data with [info panels](../../discover/info-panels.md). Datagrok's info panels contain molecule-specific functions and update dynamically and interactively.

![Toxicity, Gasteiger Partial Charges, Solubility Prediction](../../uploads/gifs/chem-model-augment.gif "Toxicity, Gasteiger Partial Charges, Solubility Prediction")

By opening, info panels automatically show the general functions:

* Structure 2D – visualizes a molecule in 2D.
* Structure 3D – visualizes a molecule in 3D through generating a MOL file.
* Molfile – provides a MOL file.

All the rest functions fall into Chemistry, Biology, Databases groups:

```mdx-code-block
<Tabs>
<TabItem value="chemistry" label="Chemistry" default>
```

Chemistry group helps explore chemical features and includes:

* [Gasteiger Partial Charges](functions/gasteiger-charges.md) - visualizes and highlights partial charges in a molecule.
* [Chem descriptors](descriptors.md) – calculates and displays the specified descriptors for a molecule.
* Properties – yields the list of physical and chemical properties: empirical formula, molecular weight, hydrogen bond acceptor (HBA) and donor (HDA) values, LogP and LogS, polar surface area (PSA), number of rotatable bonds and stereocenters, also IUPAC name.

```mdx-code-block
</TabItem>
<TabItem value="biology" label="Biology">
```

Biology group helps explore drug design related features and includes:

* [Toxicity](info-panels/toxicity-risks.md) –  predicts the toxicity scores. Consists of such categories as mutagenicity, tumorogenicity, irritating effects, reproductive effects.
* [Structural alerts](info-panels/structural-alerts.md) – highlights the fragments in structures that might greatly increase the toxicity and other problematic structural features.
* [Drug likeness](info-panels/structural-alerts.md) – predicts a score that shows how likely this molecule is to be a drug.
* ADME/Tox - performs ADME and Toxicity prediction.

```mdx-code-block
</TabItem>
<TabItem value="Databases" label="Databases" default>
```

Databases group helps access ChEmbl, PubChem and DrugBank public databases and includes:

* Identifiers - fetches all known identifiers for the specified structure across [UniChem databases](https://www.ebi.ac.uk/unichem/).
* Substructure and similarity search - provides search across ChEmbl, PubChem and DrugBank public databases.

```mdx-code-block
</TabItem>
</Tabs>
```

:::info
In addition to predefined info panels, users also can develop their own using any scripting language supported by the Datagrok platform.
:::

## Calculate

For chemical structures, you can run two types of calculation functions: descriptors/fingerprints and mapping functions.

### Descriptors and fingerprints

To vectorize molecular graph data, the platform offers [descriptors](descriptors.md) and [fingerprints](fingerprints.md).

While descriptors are more physical and fingerprints are more abstract vectors, both enable the following procedures:

* Similarity and diversity search
* Chem space dimensionality reduction
* SAR analysis
* Machine learning predictive modeling

We support the following descriptors: Lipinski, Crippen, EState, Graph, MolSurf, and [others](descriptors.md).

![descriptors](descriptors.gif "Calculate descriptors")

We support the following fingerprints: RDKFingerprint, MACCSKeys, AtomPair, TopologicalTorsion, Morgan/Circular, and [others](fingerprints.md).

<!-- GIF with fingerprints. Doesn't work yet! -->

### Mapping functions

Datagrok's mapping functions calculate molecules' textual identifiers known as [International Chemical Identifiers](https://en.wikipedia.org/wiki/International_Chemical_Identifier), interlinking the same molecules between multiple databases.

With just a few clicks, you can convert your structures to InChI and InChI keys, its hashed version.

![Mapping functions](mapping-functions.gif "Conversions to InChI and InChI keys")

## Transform

Datagrok lets transform your molecular structures for standardization and/or augmentation purposes.  
In addition to general [data wrangling](../../transform/data-wrangling.md) procedures, you can use chemistry-specific ones: curation and mutation.

### Curation

[Сhemical structure curation](chem-curate.md) standardizes your chemical structures and helps avoid data-associated errors such as duplicated vectors in the training set or incorrect structure representation.

:::tip
By chemical structure curation you can improve your SAR analysis or model's prediction accuracy.
:::

We offer the following curation methods:

* Kekulization
* Reionization
* Normalization
* Neutralization
* Tautomerization
* Main fragment selection

![Chemical structure curation](chem_curation_demo.gif)

### Mutation

Mutation function generates new structures based on the specified one. Generation is combinatorial and has the following parameters:

* `step` - number of permutations applied to the structure
* `max random result` - number of output structures

## Search

Оur search functions help navigate across molecular structures and include the following types:

* ### [Similarity and diversity search](similarity-search.md)

  For the reference structure, finds the 10 most similar/diverse ones in your dataset.  
  Similarity search sorts all the found structures by similarity score. You can try multiple [distance metrics](similarity-spe#available-distance-metrics) for different similarity scoring.

  ![Similarity search](../../uploads/gifs/similarity-search.gif "Similarity search")

* ### [Substructure search](substructure-search.md)

  Sketcher-supported function that searches for query structural pattern in the data source.

  In Datagrok, substructure search implements several advanced features:

  * Sketcher-oriented structural adjustment of the found molecules
  * Explicit hydrogen search support
  * Extended toolset for aromaticity search
  * Multicolumn filtering

:::info
In addition to an uploaded dataset, similarity and substructure search also support public databases (Chembl, PubChem, DrugBank) and [relational databases](db-substructure-similarity-search.md).
:::
<!-- #GIF с поиском субструктур в опен-сорсных БД, плохо работают -->

## Analyze

For chemical analysis, Datagrok offers multiple methods. Functionally, they fall into 2 categories:

```mdx-code-block
<Tabs>
<TabItem value="Analyze SAR" label="Analyze SAR" default>
```

Identifies structure-based bioactivity values and includes:

* [R-group analysis](r-group-analysis.md)

  Identifies R-group branches connected to a scaffold. For comparative visualizations, trellis plot is available.

  <!-- new GIF, doesn't work -->

* Activity Cliffs

  Finds similar compounds with different activity in your dataset.  
  Choose similarity percent of neighbors, run the function and explore activity cliffs displayed on 2D plots.

<!-- GIF with Activity cliffs, doesn't work -->

* Structural Alerts

  Flags potential chemical hazards taken from several rule sets containing
[1,251 substructures in the SMARTS format](https://raw.githubusercontent.com/PatWalters/rd_filters/5f70235b387baa39669f25d95079e5dfec49a47c/rd_filters/data/alert_collection.csv).
  The presence of any of these substructures triggers a structural alert shown on the info panel.

<!-- GIF with Structural Alerts, doesn't work -->

```mdx-code-block
</TabItem>
<TabItem value="Analyze Structure" label="Analyze Structure">
```

Investigates exclusively structural pecularities and includes:

* Chem space

  Explores structural similarity of your molecular data on 2D clusters. For dimensionality reduction, you can choose either t-SNE or UMAP algorithm with any of Cosine, Sokal, Tanimoto, Assimetric distance metrics.

  <!-- GIF, but tooltip bagged! -->

* Elemental Analysis
  
  Represents molecules as the sum of its constituent atoms with optional radar plots.

  ![Elemental Analysis](elemental-analysis.gif "Elemental analysis")

  :::tip

  Use elemental analysis to filter the molecules on appropriate atoms  

  ![Elemental filtering](elemental-filtering.gif "Elemental filtering")
  :::

* Scaffold Tree

  Constructs tree hierarchy based on molecular scaffolds in your dataset. Use the hierarchy to filter or select the corresponding rows in the dataset.

  ![Scaffold Tree](scaffold-tree-filter-select.gif "Scaffold Tree")

```mdx-code-block
</TabItem>
</Tabs>
```

## Predict

Datagrok   enables   machine   learning   predictive   models   by   using
chemical properties, descriptors, and fingerprints as features, and the
observed properties as targets. Within Datagrok, you can train, assess, execute,
share and use your models in pipelines.

We offer several models: PCA, Naive Bayes, K-Means, Distributed Random Forest, Gradient Boosting, XGBoost, Generalized Linear Modeling, AutoML and Deep Learning.

![Training](../../uploads/gifs/chem-train-model.gif "Training")