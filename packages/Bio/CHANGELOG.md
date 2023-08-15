# Bio changelog

## 2.9.0 (WIP)

*Dependency: datagrok-api >= 1.13.3*

### Features

* Add WebLogo property showPositionLabels

### Bug fixes 

* Fix vdRegionsViewer viewer package function name consistency
* GROK-13310: Bio | Tools: Fix Split to monomers

## 2.8.2 (2023-08-01)

This release focuses on improving the monomer cell renderer.

*Dependency: datagrok-api >= 1.13.3*

### Features

* Added sample datasets for natural and synthetic peptide sequences
* Added sample dataset for cyclic sequences with HELM notation

### Bug fixes

* GROK-13659: Bio | Tools: Fix MaxMonomerLength for Macromolecule cell renderer

## 2.8.1 (2023-07-24)

This release focuses on improving the monomer cell renderer.

*Dependency: datagrok-api >= 1.13.3*

### Features

* Monomer cell renderer now defaults to 6 characters per monomer.

## 2.8.0 (2023-07-21)

This release focuses on improving feature stability and usability.

*Dependency: datgarok-api >= 1.13.3*

### Features

* Add Copy group to Macromolecule cell context menu.
* Add DefaultSeparator package property settings.

### Bug Fixes

* Fix VdRegionsViewer filter source checkbox tooltip, workaround.

## 2.7.2 (2023-07-21)

This release focuses on improvements and bug fixes.

*Dependency: datagarok-api >= 1.13.3*

### Features

* Set default values in all dialogs where appropriate.
* Detected Helm monomer type for separator data and made it usable for MSA.
* Added alignment options to **Kalign**.
* Added separator support for **Sequence Space** and **Activity Cliffs**.
* Tooltip: shows monomer atomic structure for macromolecules.
* For macromolecule cells, added the ability to show composition ratios in the property panel.
* **Top menu**: organized the items into groups **SAR**, **Structure**, **Atomic level**, and **Search**. 

### Bug Fixes

* GROK-13048: Activity cliffs identification for macromolecules.
