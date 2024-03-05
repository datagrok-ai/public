# compute-utils changelog

## 1.25.2 (2024-03-05)

### Features

- PLV: Added additional check on step enabling
- ModelCatalog: Added help switching

## 1.25.1 (2024-03-04)

### Bug fixes

- Sensitivity Analysis: fixed input form styles

## 1.25.0 (2024-03-04)

### Features

- RFV: Now supports readme files via `meta.readme` tag
- Wizard: Now support help panel state saving 

## 1.24.1 (2024-03-02)

### Bug fixes

- ModelCatalog: Fixed error with mandatory groups

## 1.24.0 (2024-03-01)

### Features

- RFV: Customizable data upload feature

## 1.23.0 (2024-02-28)

### Features

- SA: Updated viewers
- SA: Added analysis by value in custom column

### Bug fixes

- SA: Fixed bug last row analysis

## 1.13.14 (2023-08-30)

### Features

- HistoryInput allow skip df load

## 1.13.13 (2023-08-29)

### Features

- HistoryInput redesigned

## 1.13.9 (2023-08-23)

### Fixed bugs

- Fixed bug with empty captions in RFV blocks

## 1.13.8 (2023-08-18)

### Features

- RFV: Inputs are now disabled during computations
- RFV: editState is now properly saved
- FileInput: Moved icons into options container
- Default users and groups IDs are now taken from API

## 1.13.6 (2023-08-15)

### Features

- FileInput now extends DG.InputBase
- RFV: Historical runs become non-historical on any input change
- RFV: Ribbon panel now reacts on "historicity" of the run
- RFV: Inputs now have a delay before calculation re-run
- RFV: DataFrame export to Excel now uses built-in Excel tables

### Fixed bugs

- Removed workarounds for GROK-13335 & GROK-13337
- PLV: Fixed bug with extra tabs enabled on historical run load

## 1.13.5 (2023-08-10)

### Features

- Heavily optimized export functions for RFV
- RFV now supports format setting for float inputs using `format` tag
- RFV now supports precision setting for float outputs using `precision` tag`

### Fixed bugs

- RFV now saves selected output tab on function rerun
- RFV now fully hides output during function run panel if no input tabs is provided
