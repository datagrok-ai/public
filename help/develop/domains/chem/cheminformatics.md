---
title: "Cheminformatics"
---

In addition to being a general-purpose extensible platform for scientific computing, Datagrok provides multiple options
for developing cheminformatics solutions on top of it. Depending on your needs, use one or more of the following ones,
or come up with your own solution.

* For **Typical tasks** such as calculating descriptors, substructure or similarity search, finding MCS or performing an
  R-group analysis, consider using the [Datagrok JS API](#datagrok-js-api).
* For **custom client-side computations**, consider using either [OpenChemLib.JS](#openchemlibjs),
  or [RDKit built for WebAssembly](#rdkit-in-webassembly).
* For **custom server-side computations**, a popular option is using
  [RDKit in Python](#rdkit-in-python). Python scripts can be seamlessly embedded into Datagrok
  via [Scripting](../../../compute/scripting.md).

## Datagrok JS API

For convenience, some of the commonly used functions are exposed via the [JS API](../packages/js-api.md). Use `grok.chem`
instance to invoke the methods below. Note that most of them are asynchronous, since behind the scenes they either use a
server-side assist, or a browser's separate Web Worker to call [RDKit JS](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h) functions.

```
searchSubstructure(column, pattern, settings = { substructLibrary: true });
getSimilarities(column, molecule, settings = { /* reserved */ })
findSimilar(column, molecule, settings = { limit: ..., cutoff: ... })
diversitySearch(column, metric = METRIC_TANIMOTO, limit = 10);
rGroup(table, column, core);
mcs(column);
descriptors(table, column, descriptors);
```

The three first functions, `searchSubstructure`, `getSimilarities` and `findSimilar`, are currently made client-based,
using [RDKit JS](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h). We are planning to support both server-side and browser-side
modes for these functions in the future, while the API remains as described below.

### Substructure search

`searchSubstructure(column, pattern = null, settings = { substructLibrary: true })`

This function performs substructure search using RDKit JS. It returns a [BitSet](../packages/js-api.md#bitset).
If i-th element of the input `column` contains a substructure given by a `pattern` the i-th bit set to 1, otherwise the bit is set to 0.

`column` is a column of type String. It contains molecules in any notation [supported by RDKit JS](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h): smiles, cxsmiles, molblock, v3Kmolblock, and inchi. Same applies to a `pattern`: this is a string representation of a molecule in any of the previous notations.

The `settings` object allows passing the following parameters:

* `molBlockFailover`: smarts which used as a substructure if `pattern` is invalid

To make search maximum efficient substructure search function maintains a cache of
pre-computed Pattern fingerprints (a topological fingerprint optimized for substructure screening) and specific RDKit *Mol objects*. Fingerprints are used for preliminary filtration while *Mol objects* are used for graph-based search method [`get_substruct_match`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html). Fingerprints are pre-calculated for the whole column while *Mol objects* are cached only for defined number of first molecules. One the one hand it allows to decrease memory consumption and on the other hand speeds up the search. Once the function `searchSubstructure` meets the `column` previously unmet, it builds the cache of fingerprints and *Mol objects* for the molecules in `column`. Later it is used for substructure search as long as the function is invoked for the same `column`.

Sometimes it happens that the molecule strings are invalid or not supported by RDKit JS. These inputs are skipped when indexing. For such unsupported molecule, the corresponding bit in the bitset remains as 0. Thus, the bitset is *always* of the
input `column`'s length.

### Similarity scoring

Similarity scoring functions
use [Morgan fingerprints](https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints)
and [Tanimoto similarity](https://en.wikipedia.org/wiki/Chemical_similarity) (also known as Jaccard similarity) to rank molecules from a `column` by their similarity to a given `molecule`.

Identically to the convention of `searchSubstructure`, a `column` is a column of type String containing molecules, each
in any notation RDKit JS [supports](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h): smiles,
cxsmiles, molblock, v3Kmolblock, and inchi. Same applies to a `molecule`: this is a string representation of a molecule
in any of the previous notations.

#### Find most similar molecules sorted by similarity

`findSimilar(column, molecule, settings = { limit: ..., cutoff: ... })`

The default settings are `{ limit: Number.MAX_VALUE, cutoff: 0.0 }`, thus the function ranks and sorts all molecules by
similarity.

Produces a Datagrok [`DataFrame`](https://datagrok.ai/js-api/classes/dg.DataFrame) of three
columns ([code sample](https://public.datagrok.ai/js/samples/domains/chem/similarity-scoring-sorted)):

* The 1-st column, named `molecule`, contains the original molecules from the input `column`
* The 2-nd column, named `score`, contains the corresponding similarity scores of the range from 0.0 to 1.0. The
  DataFrame is sorted descending by this column
* The 3-rd column, named `index`, contains indices of the molecules in the original input `column`

#### Scoring each molecule against a given molecule

`getSimilarities(column, molecule = null, settings = { sorted: false });`

Produces a Datagrok [`DataFrame`](https://datagrok.ai/js-api/classes/dg.DataFrame) with a
single [`Column`](https://datagrok.ai/js-api/classes/dg.Column), where the i-th element contains a similarity score for
the i-th element of the
input `Column` ([code sample](https://public.datagrok.ai/js/samples/domains/chem/similarity-scoring-scores)).

Similarly to the `searchSubstructure` function, these functions maintain a cache of
pre-computed fingerprints. Once the function `findSimilar` or `getSimilarities` meets the `column` previously unmet, it
builds the cache of computed fingerprints for the molecules in `column`. Later it is used to compute similarity
scores as long as the function is invoked for the same `column`.

To only build (aka prime) this cache without performing an actual search, e.g. in case of pre-computing the fingerprints
in the UI, call the function passing `molecule` as an empty string `''`.

In case molecule strings are not supported by RDKit JS they are skipped when indexing. For
such unsupported molecules, null is returned instead of a fingerprint. Thus, the output Column (or
DataFrame) is *always* of the input `column`'s length.
.

Examples (see on [GitHub](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples/scripts/domains/chem))

* [Calculating descriptors](https://public.datagrok.ai/js/samples/domains/chem/descriptors)
* [Substructure search](https://public.datagrok.ai/js/samples/domains/chem/substructure-search)
* [Diversity search](https://public.datagrok.ai/js/samples/domains/chem/diversity-search)
* [Most common substructure](https://public.datagrok.ai/js/samples/domains/chem/mcs)
* [R-Group analysis](https://public.datagrok.ai/js/samples/domains/chem/r-group)

### Molecule rendering

#### SVG rendering with OpenChemLib

Use `grok.chem.svgMol` to render a molecule into a `div` as an SVG, using [OpenChemLib](https://github.com/cheminfo/openchemlib-js) library:

```
ui.dialog()
  .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', 300, 200))
  .show();
```

#### Rendering to canvas with RdKit and scaffolds

Use `grok.chem.canvasMol` to render a molecule, aligned to bonds, using RdKit. To highlight a scaffold
specify it as [SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification). 
The molecule string can be specified in any format supported by RdKit. Here is a complete example:

```
let root = ui.div();
let canvas = document.createElement('canvas');
canvas.width = 300; canvas.height = 200;
root.appendChild(canvas);
await grok.chem.canvasMol(
  0, 0, canvas.width, canvas.height, canvas,
  'COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3',
  'c1ccccc1');
ui
.dialog({title:'Molecule'})
.add(root)
.show();
```

The method is currently asynchronous due to technical limitations. It will be made synchronous in the future.

Run the above examples live here: [public.datagrok.ai](https://public.datagrok.ai).

## Openchemlib.js

[OpenChemLib.JS](https://github.com/cheminfo/openchemlib-js) is a JavaScript port of the
[OpenChemLib](https://github.com/actelion/openchemlib) Java library. Datagrok currently uses it for some of the
cheminformatics-related routines that we run in the browser, such as rendering of the molecules.

Here is
an [example of manipulating atoms in a molecule](https://public.datagrok.ai/js/samples/domains/chem/mol-atoms-bonds)
using openchemlib.js.

## Rdkit in python

[RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html) are Python wrappers for RDKit, one of the best
open-source toolkits for cheminformatics. While Python scripts get executed on a server, they can be seamlessly embedded
into Datagrok via [Scripting](../../../compute/scripting.md).

Here are some RDKit in
Python-based [cheminformatics-related scripts](https://github.com/datagrok-ai/public/tree/master/packages/Chem/scripts)
in the public repository.

## Rdkit in WebAssembly

Recently, Greg Landrum, the author of RDKit, has introduced
[a way to compile its C++ code to WebAssembly](https://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html)
, thus allowing to combine the performance and efficiency of the carefully crafted C++ codebase with the ability to run
it in the browser. This approach fits perfectly with Datagrok's philosophy of performing as much computations on the
client as possible, so naturally we've jumped on that opportunity!

See also:

* [Cheminformatics](../../../datagrok/solutions/domains/chem/chem.md)
* [JavaScript development](../../develop.md)
