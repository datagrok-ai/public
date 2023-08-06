# bio changelog

## 5.35.0 (2023-08-06)

Optimize with splitterAsFastaSimple returning ISeqSplitted allowing to speed up processing of
single character alphabets.

## Breaking changes

* SplitterFunc now return type ISeqSplitted instead of string[]

### Features

* Add WebLogo property showPositionLabels

## 5.34.1 (2023-08-01)

Patch release for a small fix for Macromolecule cell renderer

*Dependency: datgarok-api >= 1.10.2*

### Bug fixes

* GROK-13659: Bio | Tools: Fix MaxMonomerLength Macromolecule cell renderer

## 5.34.0 (2023-07-21)

This release focuses on improving usability.

*Dependency: datgarok-api >= 1.10.2*

### Features

* Added `NotationConverter.getConverter()` getting a function to convert a single value.

### Breaking change

* Removed `NotationConverter.convertStringToHelm()`.
