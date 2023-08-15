# statistics changelog

## 1.2.3(2023-08-12)

### Features

* [#2105](https://github.com/datagrok-ai/public/issues/2105): Removed point and fit line color default values
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added title parameter to IFitChartOptions

## 1.2.2 (2023-08-04)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Changed minBound and maxBound to min and max in parameterBounds

## 1.2.1 (2023-08-02)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added points log in fit
  * Added droplines calculation for IC50
  * Don't render confidence intervals by default
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Wrote TS docs for the code

## 1.2.0 (2023-07-21)

This release focuses on adding new functionality and improving stability.

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

### Features

* Exposed AUC and R2 functions.
