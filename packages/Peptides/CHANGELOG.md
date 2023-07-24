# Peptides changelog

## 1.10.0 (2023-07-19)

This release focuses on improving analysis stability and usability.

*Dependency: datgarok-api >= 1.15.4*

### Features

* Improved activity distribution plots to show selected vs. all
* Added position number in main table view column headers
* The Monomer-Position viewer in the Mutation Cliffs mode now shows the number of unique sequences that have monomer mutation at the selected position
* The Monomer-Position viewer in the Mutation Cliffs mode now selects unique sequences that have monomer mutation at the selected position
* The Mutation Cliffs panel now shows unique sequences and corresponding values from columns visible in the table view

### Bug Fixes

* Fixed bug when filter couldn't be applied to the original table after analysis start
