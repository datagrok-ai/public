# Peptides changelog

## 1.10.1 (2023-07-24)

This release focuses on improving analysis stability and usability.

*Dependency: datgarok-api >= 1.15.4*

### Features

* Invariant Map is now selecting sequences instead of filtering.
* Minor UI/UX adjustments.

### Bug Fixes

* Fixed mutation pair selection resulted in wrong unique sequences in Mutation Cliffs bug.
* Fixed position column labels not showing if the page was scaled.

## 1.10.0 (2023-07-19)

This release focuses on improving analysis stability and usability.

*Dependency: datgarok-api >= 1.15.4*

### Features

* Improved activity distribution plots to show selected vs. all.
* Added position number in main table view column headers.
* The Monomer-Position viewer in the Mutation Cliffs mode now shows the number of unique sequences that have monomer mutation at the selected position.
* The Monomer-Position viewer in the Mutation Cliffs mode now selects unique sequences that have monomer mutation at the selected position.
* The Mutation Cliffs panel now shows unique sequences and corresponding values from columns visible in the table view.

### Bug Fixes

* Fixed bug when filter couldn't be applied to the original table after analysis start.
