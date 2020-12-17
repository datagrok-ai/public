<!-- TITLE: Cheminformatics Development -->
<!-- SUBTITLE: -->
 
# Cheminformatics Development

As a general-purpose extensible platform for scientific computing, Datagrok provides multiple options 
for developing custom cheminformatics-related functionality. Depending on your needs, use one or more
of the following ones, or come up with your own solution.

For performing one of the typical tasks such as calculating descriptors, substructure or similarity search,
finding MCS or performing an R-group analysis, consider using [Grok JS API](#grok-js-api).

For custom client-side computations, consider using either [OpenChemLib.JS](#openchemlib.js), or
[RDKit built for WebAssembly](#rdkit-in-webassembly).

For custom server-side computations, a popular option is using 
[RDKit in Python](#rdkit-in-python). Python scripts can be seamlessly embedded into Datagrok
via [Scripting](scripting.md).    

## Grok JS API

For convenience, some of the commonly used functions are exposed via the [JS API](js-api.md). Use `grok.chem` instance to invoke the methods below. Note that most of them are asynchronous, since behind the scenes they either use a server-side assist, or a browser's separate Web Worker to call [RDKit-JS](https://github.com/rdkit/RDKitjs) functions.

```
substructureSearch(column, pattern, settings = { substructLibrary: true });
similarityScoring(column, molecule, settings = { sorted: false });
diversitySearch(column, metric = METRIC_TANIMOTO, limit = 10);
rGroup(table, column, core);
mcs(column);
descriptors(table, column, descriptors);
```

The two first functions, `substructureSearch` and `similarityScoring`, are currently made client-based, using [RDKit-JS](https://github.com/rdkit/RDKitjs) library. We are planning to support both server-side and browser-side modes for these functions in the future, while the API remains as described below.

### Substructure Search

`substructureSearch(column, pattern, settings = { substructLibrary: true })`

The `settings` object allows passing the following parameters:

* `substructLibrary`: a boolean indicating whether to use the 'naive' search method, using [`get_substruct_match`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html) method of RDKit (`substructLibrary: false`), or (`substructLibrary`) from RDKit. The default is `true`.

`column` should be a string column containing molecules in any notation RDKit-JS [supports](https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.h): smiles, cxsmiles, molblock, v3Kmolblock, and inchi. Same applies to a `pattern`: this can be a string representation of a molecule in any format of the above.

For the chosen `substructLibrary` method, the first call to the function with a `column` previouly unseen by the `substructureSearch` function shall invoke building a library. The library is used for conducting further searches in the same given column while the `pattern` is varying.

The function's user doesn't need to manage the library, as it is done by the function automatically. As long as the function detects it has received a new column instead of the one previously seen, it will rebuild the library for this newly seen column.

It is possible to call `substructureSearch` solely for initializing the library, without an actual search. For this, call it without passing `pattern`, or passing it with the value on `null`.

Sometimes it happens that the molecule strings aren't supported by RDKit-JS. These inputs are skipped when indexing, but an error log entry is added for each such entry to the browser's console together with the actual molecule string.

### Similarity Scoring

TBW

Examples (see on [github](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples/scripts/domains/chem))
* [Calculating descriptors](https://public.datagrok.ai/js/samples/domains/chem/descriptors)
* [Substructure search](https://public.datagrok.ai/js/samples/domains/chem/substructure-search)
* [Similarity search](https://public.datagrok.ai/js/samples/domains/chem/similarity-search)
* [Diversity search](https://public.datagrok.ai/js/samples/domains/chem/diversity-search)
* [Maximum common substructure](https://public.datagrok.ai/js/samples/domains/chem/mcs)
* [R-Group analysis](https://public.datagrok.ai/js/samples/domains/chem/r-group)

## OpenChemLib.JS

[OpenChemLib.JS](https://github.com/cheminfo/openchemlib-js) is a JavaScript port of the 
[OpenChemLib](https://github.com/actelion/openchemlib) Java library. Datagrok currently uses it
for some of the cheminformatics-related routines that we run in the browser, such as
rendering of the molecules, and performing in-memory substructure search. 

Here is an [example of manipulating atoms in a molecule](https://public.datagrok.ai/js/samples/domains/chem/mol-atoms-bonds) using openchemlib.js.

## RDKit in Python

[RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html) are Python wrappers for 
RDKit, one of the best open-source toolkits for cheminformatics. While Python scripts get executed
on a server, they can be seamlessly embedded into Datagrok via [Scripting](scripting.md).  

Here are some RDKit in Python-based [cheminformatics-related scripts](https://github.com/datagrok-ai/public/tree/master/packages/ChemScripts/scripts/python)
in the public repository.

## RDKit in WebAssembly

Recently, Greg Landrum, the author of RDKit, has introduced 
[a way to compile its C++ code to WebAssembly](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html),
thus allowing to combine the performance and efficiency of the carefully crafted C++ codebase with the ability to
run it in the browser. This approach fits perfectly with Datagrok's philosophy of performing as much 
computations on the client as possible, so naturally we've jumped on that opportunity!

Here is an [example of a package](https://github.com/datagrok-ai/public/tree/master/packages/RDKitDemo) 
that exposes two functions. One of them, "rdKitInfoPanel", demonstrates how functions can be dynamically
discovered and used as an info panels. Below is the result, the "RDKitInfo" panel on the right
contains the structure that gets rendered as SVG using the RDKit library.

![](rdkit-info-panel.png)  

See also:
* [Cheminformatics](../domains/chem/cheminformatics.md)
* [JavaScript Development](develop.md)
