# ml changelog

## 6.4.13 (2024-03-07)

Add Connectivity to MCL, Harmonized macromolecule distance functions.

## 6.4.12 (2024-03-06)

### Bug fixes

* Fix `values` type to ArrayLike for DistanceMatrixService.calc

## 6.4.11 (2024-02-27)

Add support for post-processing function in dimensionality reduction. The method `Dimensionality Reduction`
(see Top Menu > ML > Dimensionality Reduction) now allows to apply post-processing functions to the resulting
embeddings.

## 6.3.39 (2023-07-21)

This release focuses on improving usability.

*Dependency: datagarok-api >= 1.15.0*

### Features

Moved the `DistanceMatrix` class from [bio](https://github.com/datagrok-ai/public/tree/master/libraries/bio) to ml lib.
