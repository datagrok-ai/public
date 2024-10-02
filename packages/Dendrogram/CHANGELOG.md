# Dendrogram changelog

## 1.2.33 (2024-09-02)

Fix for dependency math lib version update

## 1.2.32 (2024-08-23)

### Bug fixes

* GROK-15949:Fix open GridWithTree: error
* Fix error on second Hierarchical Clustering
* Fix error on mouse hover and current row changed

## 1.2.31 (2024-08-19)

Fix Hierarchical clustering namings.

## 1.2.30 (2024-07-23)

Bump dependencies versions, datagrok-api to 1.20.0

## 1.2.29 (2024-04-23)

Bump dependencies versions, datagrok-api to 1.18.0, gridext to 1.3.71
GROK-15153: Fix test Demo.heatMapDemo

## 1.2.28 (2024-03-30)

### Features

* #2707: Add original and canonical to monomer

## 1.2.27 (2024-03-07)

### Bug fixes

* Bump dependencies version, NEWICK_EMPTY

## 1.2.26 (2024-03-06)

### Bug fixes

* Bump dependencies version
* Fix tree traverse order, fix tests

## 1.2.18 (2023-07-21)

This release focuses on improving analysis stability and usability.

* Dependency: datagarok-api >= 1.15.0*

### Features

* Implemented the ability to select all leaves from a certain node.
* Added a separate loader view to the dendrogram.
* Added the ability to reset zoom for Dendrogram.
* Moved hierarchical clustering script to the client side.
* Allowed to switch distance calculation method for macromolecules (Hamming or Levenstein) in case of MSA.
* Enabled Dendrogram to work with semType `macromolecules`.
