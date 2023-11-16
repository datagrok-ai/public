# bio changelog

## 5.39.0 (2023-10-25)

### Features

* Add VdRegionsViewer `filterSource` property.
* Add WebLogo `valueAggrType`, `valueColumnName` properties.
* Add routines to get/set user monomer library settings.
* Add Molecule3DUnitsHandler

### Bug fixes

* Fix getSplitterWithSeparator for empty seq.
* Fix MonomerPlacer getPosition for empty seq.
* Move maxMonomerLength from column's temp to tags
* Add property fitWidth for VdRegionsViewer
* Fix UnitsHandler to not replace existing tags
* Fix VdRegionsViewer fit width accounting position margin of WebLogo
* Fix UnitsHandler to allow empty alphabet if annotated with tag.
* Fix UnitsHandler.maxLength to be zero for empty data.
* Fix monomer lib handler for PolyTool
* Fix splitter with separator for quoted gaps
* Fix default value for pdbTag prop of Biostructure

## 5.38.0 (2023-09-06)

### Features

* Add UnitsHandler `maxLength`, `posList`, `splitter` properties.
* Add UnitsHandler `getJoiner`, `getConverter` methods, granulate joiners and converters methods.

### Breaking changes

* Remove NotationConverter class, merge functionality into UnitsHandler.

## 5.37.0 (2023-08-30)

### Features

* Added `calculateSimilarity` and `calculateIdentity` functions.
* Added `calculateScores` function.

## 5.36.1 (2023-08-18)

### Bug fixes

* Restore utils `getSplitter` method, it is required for MacromoleculeDifference column.

## 5.36.0 (2023-08-10)

This release introduces chemical sequence similarity functionality.

### Features

* Added `sequenceChemSimilarity` function.

## 5.35.0 (2023-08-06)

Optimize with splitterAsFastaSimple returning ISeqSplitted allowing to speed up processing of
single character alphabets.

## Breaking changes

* SplitterFunc now returns `ISeqSplitted` instead of `string[]`.

### Features

* Add WebLogo property showPositionLabels.

## 5.34.1 (2023-08-01)

Patch release for a small fix for the Macromolecule cell renderer.

### Bug fixes

* GROK-13659: Bio | Tools: Fix MaxMonomerLength Macromolecule cell renderer.

## 5.34.0 (2023-07-21)

This release focuses on improving usability.

### Features

* Added `NotationConverter.getConverter()` getting a function to convert a single value.

### Breaking change

* Removed `NotationConverter.convertStringToHelm()`.
