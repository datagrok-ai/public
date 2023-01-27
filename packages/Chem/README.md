# Cheminformatics support

This package provides first-class [cheminformatics support](https://datagrok.ai/cheminformatics) for the Datagrok platform.
See it in action on [YouTube](https://www.youtube.com/watch?v=k1NVdTRpYOM&ab_channel=Datagrok).

The Datagrok platform provides rich exploratory data analysis capabilities, advanced data mining
techniques such as multivariate analysis and stochastic proximity embedding, and out-of-the-box support for
predictive modeling and scientific computations. By providing first-class 
[cheminformatics support](https://datagrok.ai/cheminformatics), this package turns Datagrok into 
a comprehensive platform for working with chemical and biological data.
 
Regarding the performance, our goal is to be able to open chemical datasets of up to 10 millions small molecules completely 
in the browser, and interactively perform commonly used operations such as substructure and similarity search 
without having to rely on a server. In order to hit these goals, we are using a couple of techniques. First of all, 
we are leveraging Datagrok's capability to [efficiently work with relational data](https://datagrok.ai/help/develop/performance).
For cheminformatics, we are relying on the [RDKit library compiled to WebAssembly](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html). 
This not only gives us the ability to execute C++ code at the native speed, but also take full advantage of the 
modern multicore CPUs by running computations in multiple threads. 

Here are some of the features:

* Molecule sketching
  * OpenChemLib
* SAR analysis
* Chembl integration
* Pubchem integration
* Works completely on the client side where possible.
* [RDKit](https://www.rdkit.org)-based rendering. 
  * Highlighting substructures on search (in progress) 
  * Aligning to scaffolds (in progress) 
  * Rendering options (in progress)
* Fingerprint-based similarity and diversity analyses (see [video](https://www.youtube.com/watch?v=wCdzD64plEo&ab_channel=Datagrok))
* Efficient in-memory substructure and similarity searching

The following features are still in the core, but we plan to move them out to this package:

* Property calculators (server-side)
* 3D: coordinate calculation using RDKit, rendering using [NGL Viewer](http://nglviewer.org/)
* "Sketch-to-predict": run predictive models as you sketch the molecule

# Sketcher

Supports multiple sketchers (MarvinJS, OpenChemLib, Ketcher).

## Favorite and recent structures

Access the recently sketched structures from the `☰ -> Recent` menu.

`☰ -> Favorites` contains your favorite structures. To add current molecule
to the favorites, click on `☰ -> Favorites -> Add to Favorites`.

## Molecule queries

Out-of-the-box, you can paste SMILES, MOLBLOCK, and InChi keys into the input field, and the sketcher
automatically translates it to a structure. In addition to that, you can make sketcher understand
other structure notations (such as from your company's internal database of structures) by registering
a function annotated in a special way. The following example provides support for Chembl. The important
tags are:
* `--meta.role: converter`: indicates that such a function serves as a value converter
* `--meta.inputRegexp: (CHEMBL[0-9]+)`: RegExp that is evaluated to check if this function
  is applicable to the user input. The captured group (in this case the whole input) is then
  passed to this function as a parameter. 
* `--output: string smiles { semType: Molecule }`: should return string with the semType `Molecule`

```
--name: chemblIdToSmiles
--meta.role: converter
--meta.inputRegexp: (CHEMBL[0-9]+)
--connection: Chembl
--input: string id = "CHEMBL1185"
--output: string smiles { semType: Molecule }
select canonical_smiles from compound_structures s
join molecule_dictionary d on s.molregno = d.molregno
where d.chembl_id = @id
```

![](help/molecule-queries.gif)

A molecule query does not have to be a database query, any [function](../../help/datagrok/functions/function.md)
will do. For instance, InChi query is implemented as a Python script. 

## Elemental analysis

Elemental Analysis lets you extract, visualize, and compare atom counts for the molecules in your dataset.

With this tool you can:
* get columns with atom counts (same order as in the periodic table);
* get a radar chart for every single molecule which makes it possible to compare each molecule with the whole dataset.

To run Elemental Analysis, select Chem | Elemental Analysis from the main menu.

Radar chart can be visualized in 2 ways:
1. A new column with radar charts is added.

![](help/radar-chart-grid.gif)

2. Radar chart appears in a new window. It changes automatically when you click on the new row.

![](help/radar-chart-view.gif)

## Similarity search
Toll that help scientists analyze a collection of molecules in terms of molecular similarity. It is based on applying different distance metrics (such as Tanimoto) to fingerprints.

To run Similarity search select Chem | Similarity search from the top menu.

![](help/similarity_search.png)

To change target molecule select the row with required molecule in the initial dataframe. Similarities will be recalculated.
Target molecule can be also changed by using `edit` button in the top left corner of the molecule pane.

![](help/similarity_search_sketch_target.gif)

Using property panel you can change search metrics like similarity cut off, fingerprints type or distance metric.

By selecting columns from `Molecule Properties` field you can add any fields present in your dataframe to similarity panes. Please note that if color coding is applied to a selected column it will be copied on a similarity pane. You can select to highlight background or text in the corresponding field. Ability to add fields to similarity panes simplifies analysis since multiple molecules characteristics can be easily assesed at once.

![](help/similarity_search_add_fields.gif)


See also: 
  * [Cheminformatics predictive modeling](https://datagrok.ai/help/domains/chem/chem-predictive-modeling)
  * [Datagrok JS API](https://datagrok.ai/help/develop/js-api)
  * [Datagrok Cheminformatics API](https://datagrok.ai/help/develop/cheminformatics-development)
  * [Introducing new RDKit JavaScript wrappers](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html)
  * [Webassembly](https://webassembly.org/)