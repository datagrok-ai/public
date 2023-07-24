# Bio changelog

## 2.8.0 (2023-07-21)

This release focuses on improving feature stability and usability.

*Dependency: datgarok-api >= 1.13.3*

### Features

* Add Copy group to Macromolecule cell context menu
* Add DefaultSeparator package property settings

### Bug Fixes

* Fix VdRegionsViewer filter source checkbox tooltip, workaround

## 2.7.2 (2023-07-21)

This release focuses on improving analysis stability and usability.

*Dependency: datagarok-api >= 1.13.3*

### Features

* Set default values in all dialogs where appropriate.
* Detected Helm monomer type for separator data and made it usable for MSA.
* Added alignment options to **Kalign**.
* Added separator support for **Sequence Space** and **Activity Cliffs**.
* Implemented showing monomer atomic structure in tooltips for macromolecules.
* For macromolecule cell, added the ability to show composition ratios in property panel.
* We have structured the top menu by organizing the items into groups: SAR, Structure, Atomic level, and Search. 

### Bug Fixes

* GROK-13048: Activity cliffs identification for macromolecules.
