# Bio changelog

## 2.10.0 (2023-09-06)

### Features

* GetRegion for Macromolecule.
* Top menu Bio | Convert | GetRegion dialog.
* Add tests for UnitsHandler.getRegion.
* Maintain `.positionNames` tag for GetRegion derived column.
* Use `.regions` tag annotation for GetRegion dialog.
* Fix mistyping top menu path for Identity scoring.

### Bug fixes

* Fixed UnitsHandler.posList length.

## 2.9.0 (2023-08-30)

### Features

* WebLogo: add property `showPositionLabels`.
* WebLogo: optimized with `splitterAsFastaSimple`.
* WebLogo: disable `userEditable` for `fixWidth`.
* VdRegionsViewer: optimized preventing rebuild on `positionWidth` changed and resize.
* VdRegionsViewer: to fit WebLogo enclosed on `positionWidth` of value 0.
* Introduced sequence identity and similarity scoring.

### Bug fixes

* Fix vdRegionsViewer viewer package function name consistency.
* GROK-13310: Bio | Tools: Fix Split to monomers for multiple runs.
* GROK-12675: Bio | Tools: Fix the Composition dialog error on the selection column.
* Allow characters '(', ')', ',', '-', '_' in monomer names for fasta splitter.
* WebLogo: Fix horizontal alignment to the left while `fixWidth``.
* WebLogo: Fix layout for `fixWidth`, `fitArea`, and normal modes.
* VdRegionsViewer: Fix postponed rendering for tests.
* MacromoleculeDifferenceCellRenderer: Fix to not use `UnitsHandler`.

## 2.8.2 (2023-08-01)

This release focuses on improving the monomer cell renderer.

*Dependency: datagrok-api >= 1.13.3*

### Features

* Added sample datasets for natural and synthetic peptide sequences.
* Added sample dataset for cyclic sequences with HELM notation.

### Bug fixes

* GROK-13659: Bio | Tools: Fix MaxMonomerLength for Macromolecule cell renderer.

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
