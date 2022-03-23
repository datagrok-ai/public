<!-- TITLE: Diversity search -->
<!-- SUBTITLE: -->

# Diversity search

We have implemented few tools that help scientists analyze a collection of molecules in terms of molecular similarity.
Both tools are based on applying different distance metrics
(such as Tanimoto) to fingerprints.

* [Similarity Search](similarity-search.md) - finds structures similar to the specified one
* [Diversity Search](diversity-search.md) - finds 10 most distinct molecules

These tools can be used together as a collection browser. 'Diverse structures' window shows different classes of
compounds present in the dataset; when you click on a molecule representing a class, similar molecules will be shown in
the 'Similar structures' window.

![similarity-search](similarity-search.gif)

See also:

* [Descriptors](descriptors.md)
* [Molecular fingerprints](fingerprints.md)
* JS API: [Diversity search](https://public.datagrok.ai/js/samples/domains/chem/diversity-search)
