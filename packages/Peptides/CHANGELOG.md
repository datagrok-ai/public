# Peptides changelog

## 1.12.0 (WIP)

### Features

* Improved Identity panel UI.
* Added sequence similarity functionality.
* Moved Similarity & Identity to Actions panel.
* Other UI/UX improvements.
* Peptides moved to Bio > Analyze > SAR.
* Added tooltips for all inputs and buttons.
* Added Selection panel which shows the current selection as a separate table.

### Bug Fixes

* Fixed table view tooltips.
* Fixed selection with WebLogo in the column header.
* Fixed identity scoring formula.
* Fixed table view couldn't be filtered by monomer position columns.
* Fixed identity panel not considering filtered table view.
* Fixed the Identity button not disabling on incorrect input.

## 1.11.0 (2023-08-03)

This release introduces sequence identity functionality and some stability and usability improvements.

### Features

* Added sequence identity functionality.
* UI/UX improvements.
* Added activity column choice in settings.

## 1.10.4 (2023-08-01)

### Features

* UI/UX improvements.

### Bug Fixes

* Fixed Monomer-Position viewer not responding to selection correctly.

## 1.10.3 (2023-07-31)

### Bug Fixes

* Fixed wrong values in the invariant map.
* Fixed Mutation Cliffs pair selection would still result in wrong sequences sometimes.

## 1.10.2 (2023-07-24)

### Bug Fixes

* Made invisible Mutation Cliffs pairs columns, containing sequence indexes in the original table view.

## 1.10.1 (2023-07-24)

### Features

* Invariant Map is now selecting sequences instead of filtering.
* Minor UI/UX adjustments.

### Bug Fixes

* Fixed mutation pair selection resulted in wrong unique sequences in the Mutation Cliffs bug.
* Fixed position column labels not showing if the page was scaled.

## 1.10.0 (2023-07-19)

This release focuses on improving analysis stability and usability.

### Features

* Improved activity distribution plots to show selected vs. all.
* Added position number in main table view column headers.
* The Monomer-Position viewer in the Mutation Cliffs mode now shows the number of unique sequences that have monomer mutation at the selected position.
* The Monomer-Position viewer in the Mutation Cliffs mode now selects unique sequences that have monomer mutation at the selected position.
* The Mutation Cliffs panel now shows unique sequences and corresponding values from columns visible in the table view.

### Bug Fixes

* Fixed bug when filter couldn't be applied to the original table after analysis started.
