# statistics changelog

## 1.2.0 (2023-07-21)

This release focuses on adding new functionality and improving stability.

*Dependency: datagarok-api >= 1.16.0*

### Features

* Added fit function determination.
* Added the capability to run user-defined Javascript fitting functions with caching.
* Added box plot statistics.
* Added `candlesticks` and `both` render parameter types.
* Added different marker types.
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added interface for custom parameters.
  * Added the capability to control parameter bounds.

### Bug Fixes

* Fixed inputs in the property panel.
* Fixed confidence intervals and statistics methods.
* Fixed confidence interval stroke color opacity.
* Fixed fit function creation.

## 1.1.8 (2023-05-21)

This release focuses on adding new functionality.

*Dependency: datagarok-api >= 1.13.11*

### Features

* Exposed AUC and R2 functions.
