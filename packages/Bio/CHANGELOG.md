# Bio changelog

## 2.12.3 (2024-04-03)

Updated version of openchemlib in dependencies

## 2.12.2 (2024-04-03)

Harmonized MM distance function with monomer similarity matrices.

## 2.12.1 (2024-04-02)

### Bug fixes

* Fix bioSubstructureFilter with two columns
* Fix bioSubstructureFilter error on filters reopen
* GROK-15292: Fix bioSubstructureFilter for reset

## 2.12.0 (2024-03-30)

### Features

* #2707: Add original and canonical to monomer

## 2.11.42 (2024-03-27)

### Bug fixes

* Fix MacromoleculeColumnWidget error with WebLogo disabled
* Fix WebLogo for GAP_SYMBOL
* Fix CompositionAnalysisWidget for gaps

## 2.11.41 (2024-03-26)

### Features

* Polytool: ability to use special engine to create molV3000 with CFG flags in the atoms block

## 2.11.40 (2024-03-22)

### Features

* Polytool rules file handling

## 2.11.39 (2024-03-11)

### Bug fixes

* GROK-15150: Fix display hidden/showed inputs

## 2.11.38 (2024-03-08)

### Bug fixes

* GROK-14910: PepSeA verbose output
* MSA ensures docker container for PepSeA
* Sample files harmonized

## 2.11.37 (2024-03-08)

### Bug fixes

* Check Bio publishing, PepSeA docker

## 2.11.36 (2024-02-28)

### Bug fixes

* GROK-15086: Fix GetRegion, Notation: result column does not render

## 2.11.35 (2024-02-26)

### Features

* #2706: Polytool: init rule based generation

## 2.11.34 (2024-02-20)

### Bug fixes

* Downgrade datagrok-api dependency version to 1.17.4

## 2.11.32 (2024-02-20)

### Bug fixes

* GROK-11982: Bio: Fix duplicates WebLogo on layout, test
* GROK-11983: Bio: Fix duplicates WebLogo on project, test

## 2.11.31 (2024-02-19)

### Features

* GROK-14230: Bio: Add basic UI for monomer lib files adding / validation

## 2.11.30 (2024-02-15)

### Features

* GROK-14598: Bio: Substructure filter sync between cloned views, tests

### Bug fixes

* GROK-14916: Bio: Fix biosubstructure filter for sequences of Helm
* GROK-14913: Bio: Fix To Atomic Level for sequences with gaps, tests

## 2.11.28 (2024-02-07)

### Features

### Bug fixes

* Fix detectMacromolecule allowing double quoted sequences and gaps.
* Fix for min seq length 10, tests.
* Fix user library settings for tests
* Fix test substructureFilter/helm
* Fix README.md images, GIFs size

## 2.11.0 (2023-10-25)

### Features

* Add VdRegionsViewer `filterSource` property.
* Add ToAtomicLevel for non-linear HELM structures.
* Add WebLogo aggregation function.
* Add WebLogo position tooltip with composition table (for count).
* Add PolyTool with Helm2Molfile support

### Bug fixes

* Fix GetRegion to detect semantic type and renderer for created column.
* Fix GetRegion dialog column name field for default value.
* Fix WebLogo viewer for gaps with Helm.
* Fix cell renderer for empty values MSA.
* Fix detectMacromolecule to ignore empty seqs.
* Add property fitWidth for VdRegionsViewer.
* Fix VdRegionsViewer fit width accounting position margin of WebLogo.
* Fix similarity/diversity viewer tests.
* Fix reset filters for Substructure Search filter.
* Fix Macromolecule cell renderer width limit for `devicePixelRatio` less than 1.
* Fix VdRegionsViewer `positionHeight` transmit to enclosed WebLogo.
* Fix added by Split to monomers columns to be not used as a filter.
* Fix WebLogo with `filterSource` of `Selection` to work from project.
* Fix WebLogo with `filterSource` of `Selection` display all in case of empty selection.
* Fix VdRegionsViewer for initial `filterSource`.
* Fix Macromolecule cell render to skip on closed grid.
* Fix WebLogoViewer for empty column (annotated).
* Fix WebLogoViewer and VdRegionsViewer deadlock.
* Fix Activity Cliffs error on Helm dataset
* Fix Macromolecule tooltip, context widget for detach
* Fix WebLogo optimize with postponed update positions
* Fix error in bioSubstructureFilter for Helm

## 2.10.0 (2023-09-06)

### Features

* GetRegion for Macromolecule:
  * Added dialog Top menu Bio | Convert | GetRegion.
  * Added maintaining `.positionNames` tag for GetRegion derived column.
  * Added using `.regions` tag annotation for GetRegion dialog.

### Bug fixes

* Fixed UnitsHandler.posList length.
* Fixed mistyping top menu path for Identity scoring.
* Fixed Get Region missed in short top menu.
* Fixed tests detectMacromoleculeBenchmark to specify test failed.

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
