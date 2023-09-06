# bio changelog

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
