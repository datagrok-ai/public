# Cheminformatics support (beta).

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

* Works completely on the client side where possible.
* [RDKit](https://www.rdkit.org)-based rendering. 
  * Highlighting substructures on search (in progress) 
  * Aligning to scaffolds (in progress) 
  * Rendering options (in progress)
* Fingerprint-based similarity and diversity analyses (see [video](https://www.youtube.com/watch?v=wCdzD64plEo&ab_channel=Datagrok))
* Efficient in-memory substructure and similarity searching

The following features are still in the core, but we plan to move them out to this package:

* Molecule sketching
  * OpenChemLib
* SAR analysis
* Property calculators (server-side)
* 3D: coordinate calculation using RDKit, rendering using [NGL Viewer](http://nglviewer.org/)
* Chembl integration
* Pubchem integration
* "Sketch-to-predict": run predictive models as you sketch the molecule

See also: 
  * [Cheminformatics predictive modeling](https://datagrok.ai/help/domains/chem/chem-predictive-modeling)
  * [Datagrok JS API](https://datagrok.ai/help/develop/js-api)
  * [Datagrok Cheminformatics API](https://datagrok.ai/help/develop/cheminformatics-development)
  * [Introducing new RDKit JavaScript wrappers](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html)
  * [Webassembly](https://webassembly.org/)