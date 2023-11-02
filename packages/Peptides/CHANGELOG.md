# Peptides changelog

## 1.16.0 (WIP)

### Features

* Added Select All and Desect All functionality to all viewers.
* Added Mean activity column for the Most Potent Residues viewer.
* Added Mean activity to tooltips.

### Bug Fixes

* Fixed Logo Summary Table tooltip wouldn't show if analysis table doesn't contain column names 'Cluster'.
* Fixed tooltips would show irrelevant info.

## 1.15.3 (2023-10-26)

### Bug Fixes

* Fixed docking sequence space would completely remove the viewer.

## 1.15.2 (2023-10-26)

### Features

* Improved Sequence space default parameters (changed threshold similarity from 0.8 to 0.3)

## 1.15.1 (2023-10-19)

### Features

* Color-coding Sequence space by scaled activity.
* Optimized Mutation Cliffs calculations to use workers.
* Hid WebLogo positions in Logo Summary Table.

### Bug Fixes

* Fixed Invariant Map color coding wouldn't change when changing color column.

## 1.15.0 (2023-10-12)

### Features

* Added Sequence Space viewer.
* Selection panel columns now inherti width from main table.
* Invarian Map cell values now rendered with contrast color relative to background.

### Bug Fixes

* Fixed WebLogo in header margin width.
* Widgets in Context panel will now fill all the available width.
* Fixed LST rows not resizable.
* Fixed selection inconsistency.
* Fixed Most Potent Residues keyboar navigation.
* Fixed analysis not starting when there are columns of the same name.

## 1.14.1 (2023-09-28)

### Bug Fixes

* Fixed a bug when filtering single cluster would break Logo Summary Table.

## 1.14.0 (2023-09-28)

### Features

* Added keyboard navigation to Logo Summary Table, Monomer-Position, Most Potent Residues viewers and Mutation Cliff pairs panel.
* Added icon to expand grids in property panel to fullscreen.

### Bug Fixes

* Fixed Distribution panel labels.

## 1.13.3 (2023-09-20)

### Features

* WebLogo now reacts to filtering.

### Bug Fixes

* Fixed distribution panel throwing error when selection matches filter.

## 1.13.2 (2023-09-14)

### Bug Fixes

* Downgraded `datagrok-api` dependency to version 1.16.4.

## 1.13.1 (2023-09-13)

### Bug Fixes

* Most Potent Residues: Fixed issue when the viewer failed to construct in some cases.

## 1.13.0 (2023-09-13)

### Features

* Monomer-Position and Most Potent Residues: Circle size is now based on the absolute value of mean difference.
* Scaling default colors (blue for low values and red for high) in Monomer-Position viewer.
* Highlight rows on mouse over Monomer-Position, Most Potent Residues, Logo Summary Table and WebLogo in the header.
* Changed viewer interactivity to match the platform defaults. Learn more about viewer hotkeys [here](https://datagrok.ai/help/visualize/viewers/#selection).
* Accordion with actions, mutation cliffs pairs, distribution and selection panels now shows on any selection.
* Removed steps from tutorial.
* Ignore rows with missing values when starting the analysis.

### Bug Fixes

* Logo Summary Table: Fixed interaction.

## 1.12.0 (2023-08-30)

### Features

* Added sequence identity and similairty functionality.
* Peptides moved to Bio > Analyze > SAR.
* Added tooltips for all inputs and buttons.
* Added Selection panel which shows the current selection as a separate table.
* UI/UX improvements.

### Bug Fixes

* Fixed table view tooltips.
* Fixed selection with WebLogo in the column header.
* Fixed table view couldn't be filtered by monomer position columns.

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
