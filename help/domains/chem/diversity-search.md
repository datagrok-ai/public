---
title: "Diversity search"
---

We have implemented few tools that help scientists analyze a collection of molecules in terms of molecular similarity.
Both tools are based on applying different distance metrics
(such as Tanimoto) to fingerprints.

* [Similarity Search](similarity-search.md) - finds structures similar to the specified one
* [Diversity Search](diversity-search.md) - finds 10 most distinct molecules

These tools can be used together as a collection browser. 'Diverse structures' window shows different classes of
compounds present in the dataset; when you click on a molecule representing a class, similar molecules will be shown in
the 'Similar structures' window.

To run Diversity search select Chem | Search | Diversity search from the top menu.

Use context panel to change search metrics like fingerprints type or distance metric.
Additional fields can be added to molecule panes using `Molecule Properties` field on a context panel. Color coding from
initial dataframe will be saved.

![diversity-search](img/diversity_search.png)

See also:

* [Descriptors](descriptors.md)
* [Molecular fingerprints](fingerprints.md)
* JS API: [Diversity search](https://public.datagrok.ai/js/samples/domains/chem/diversity-search)
