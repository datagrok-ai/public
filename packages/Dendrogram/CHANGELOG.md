# Dendrogram changelog

## 1.2.18 (2023-07-21)

This release focuses on improving analysis stability and usability.

*Dependency: datagarok-api >= 1.15.0*

### Features

* Implemented the ability to select all leaves from a certain node.
* Added a separate loader view to the dendrogram.
* Added the ability to reset zoom for Dendrogram.
* Moved hierarchical clustering script to the client side.
* Allowed to switch distance calculation method for macromolecules (Hamming or Levenstein) in case of MSA.
* Enabled Dendrogram to work with semType `macromolecules`.
