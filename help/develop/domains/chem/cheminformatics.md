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
  via [Scripting](../../../compute/scripting/scripting.mdx).

## Datagrok JS API

We expose some commonly used functions via [JS API](../../packages/js-api.md). To invoke the methods, use `grok.chem`. 
Most of the methods are asynchronous. Behind the scenes they either use a
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

`searchSubstructure`, `getSimilarities` and `findSimilar` are client-based functions and
use [RDKit JS](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h). We plan to support both server-side and browser-side modes. The API remains the same.

### Substructure search

`searchSubstructure(column, pattern = null, settings = { substructLibrary: true })`

This function performs substructure search using RDKit JS. It returns a [BitSet](../../packages/js-api.md#bitset).
If the i-th element in the input `column` has the pattern's substructure, the i-th bit is set to 1; otherwise, it is set to 0.

`column` stands for [column](https://datagrok.ai/js-api/classes/dg.Column) that contains molecules in any
notation [supported by RDKit JS](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h): smiles, cxsmiles, molblock, v3Kmolblock, and inchi. Same applies to a `pattern`: this is a string representation of a molecule in any of the previous notations.

The `settings` object allows passing the following parameters:

* `molBlockFailover`: smarts which used as a substructure if `pattern` is invalid

To optimize the substructure search we uses a cache of pre-computed [Pattern fingerprints](https://www.rdkit.org/docs/RDKit_Book.html#additional-information-about-the-fingerprints) and specific RDKit *Mol objects*. Fingerprints are used for preliminary filtration while *Mol objects* are used for graph-based search method [`get_substruct_match`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html). Fingerprints are calculated for each molecule in `column` and *Mol objects* are created only for defined number of first molecules. It allows to decrease memory consumption and speed up the search. When the substructure search function encounters a new `column`, it creates a cache of fingerprints and *Mol objects* for the molecules in it. Substructure search uses this cache when the function is called for the same column.

Molecule strings may be invalid or not supported by RDKit JS. During indexing we skip these inputs. The corresponding bit in the bitset remains 0.. Thus, the bitset is *always* of the input `column`'s length.

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

* `molecule` column contains the original molecules from the input `column`.
* `score` column contains the corresponding similarity scores of the range from 0.0 to 1.0. The
  DataFrame is sorted descending by this column.
* `index` column contains indices of the molecules in the original input `column`.

#### Scoring each molecule against a given molecule

`getSimilarities(column, molecule = null, settings = { sorted: false });`

Produces a Datagrok [`DataFrame`](https://datagrok.ai/js-api/classes/dg.DataFrame) with a
single [`Column`](https://datagrok.ai/js-api/classes/dg.Column), where the i-th element contains a similarity score for
the i-th element of the
input `Column` ([code sample](https://public.datagrok.ai/js/samples/domains/chem/similarity-scoring-scores)).

Similarly to the `searchSubstructure` function, these functions maintain a cache of
pre-computed fingerprints. If `findSimilar` or `getSimilarities` function encounters a new `column`, it creates a cache of fingerprints for the molecules in it. And uses this cache to compute similarity scores for the same column.

To only build (aka prime) this cache without performing an actual search, e.g. in case of pre-computing the fingerprints
in the UI, call the function passing `molecule` as an empty string `''`.

If molecules are not supported by RDKit JS, we skip them during indexing. Null is returned instead of a fingerprint. Thus, the output Column (or DataFrame) is *always* of the input `column`'s length.

Examples (see on [GitHub](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples/scripts/domains/chem))

* [Calculating descriptors](https://public.datagrok.ai/js/samples/domains/chem/descriptors)
* [Substructure search](https://public.datagrok.ai/js/samples/domains/chem/substructure-search)
* [Diversity search](https://public.datagrok.ai/js/samples/domains/chem/diversity-search)
* [Most common substructure](https://public.datagrok.ai/js/samples/domains/chem/mcs)
* [R-Group analysis](https://public.datagrok.ai/js/samples/domains/chem/r-group)

### Molecule rendering

#### SVG rendering with OpenChemLib

To render a molecule into a `div` as an SVG, use `grok.chem.svgMol`. This method uses [OpenChemLib](https://github.com/cheminfo/openchemlib-js) library.

```
ui.dialog()
  .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', 300, 200))
  .show();
```

#### Rendering to canvas with RdKit and scaffolds

To render a molecule, aligned to bonds, use `grok.chem.canvasMol`. This method uses RDKit library.
To highlight a scaffold specify it as [SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification).
You can specify the molecule string in any format supported by RdKit. Here is a complete example:

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
into Datagrok via [Scripting](../../../compute/scripting/scripting.mdx).

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
