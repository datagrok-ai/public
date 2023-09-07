# Compute changelog

## 1.14.1 (2023-09-06)

### Fixed bugs

- Fixed color and numeric values parsing in RFV

## 1.13.12 (2023-08-25)

### Fixed bugs

- Fixed rare bug with historical runs in PLV
- Fixed bug with viewers' captions duplication

## 1.13.9 (2023-08-23)

### Fixed bugs

- Fixed bug with empty captions in RFV blocks

## 1.13.8 (2023-08-18)

### Features

- RFV: Inputs are now disabled during computations
- RFV: editState is now properly saved
- FileInput: Moved icons into options container
- Default users and groups IDs are now taken from API

### Bug fixes

- History panel CSS fixed

## 1.13.6 (2023-08-15)

### Features

- FileInput now extends DG.InputBase
- RFV: Historical runs become non-historical on any input change
- RFV: Ribbon panel now reacts on "historicity" of the run
- RFV: Inputs now have a delay before calculation re-run
- RFV: DataFrame export to Excel now uses built-in Excel tables

### Bug fixes

- Removed workarounds for GROK-13335 & GROK-13337
- PLV: Fixed bug with extra tabs enabled on historical run load

## 1.13.5 (2023-08-10)

### Features

- Heavily optimized export functions for RFV
- RFV now supports format setting for float inputs using `format` tag
- RFV now supports precision setting for float outputs using `precision` tag`

### Bug fixes

- RFV now saves selected output tab on function rerun
- RFV now fully hides output during function run panel if no input tabs is provided

## 1.13.3 (2023-08-01)

- ModelCatalog now shows model readme on card click

## 1.13.2 (2023-07-28)

- RichComputationViewEditor now uses choiceInputs for string properties

## 1.13.0 (2023-07-28)

- RichComputationViewEditor now supports exporting plots' images as separate export option

## 1.12.3 (2023-07-24)

- RichComputationViewEditor now supports inputs defined in other DG functions
- Fixed package augmentation for manually created scripts
