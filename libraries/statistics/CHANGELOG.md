# statistics changelog

## 1.4.0 (2024-08-23)

### Features

* Improved fitting algorithm - didn't fit properly because was stuck in the bounds

## 1.3.6 (2024-08-13)

### Features

* Added exponential fit function

## 1.3.5 (2024-05-23)

### Bug fixes

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Changed labels (Marker Type -> Marker and Outlier Marker Type -> Outlier Marker)

## 1.3.4 (2024-05-22)

### Bug fixes

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Rollback from marker to markerType

## 1.3.3 (2024-05-22)

### Features

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Changed markerType to marker

## 1.3.2 (2024-05-21)

### Features

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Added outlier marker option

## 1.3.1 (2024-05-21)

### Features

* [#2797](https://github.com/datagrok-ai/public/issues/2797): Added log-linear fit function

## 1.3.0 (2024-05-11)

### Features

* [#2797](https://github.com/datagrok-ai/public/issues/2797): Changed optimizer for fitting

## 1.2.13 (2024-04-18)

### Features

* [#2754](https://github.com/datagrok-ai/public/issues/2754): Implemented capability just to connect the points (without fitting)
* ODEs: prepare for fitting attachment
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improve curve properties and rendering:
  * Added mergeSeries property
  * Added column labels

### Bug Fixes

* Fixed property panel names.
* [#2103](https://github.com/datagrok-ai/public/issues/2103): Renamed Show Statistics to Statistics

## 1.2.12 (2024-03-14)

### Features

* [#2105](https://github.com/datagrok-ai/public/issues/2105): Added outlier color property.

## 1.2.11 (2023-11-17)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added allowXZeroes property for logarithmic data.

## 1.2.10 (2023-11-09)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed curves not rendering because of undefined x if branch.
  * Fixed linear fit not rendering.

## 1.2.9 (2023-11-02)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Implemented linear function.

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed AUC calculation for logarithmic values.
  * Fixed exception if no x or y coordinates present.

## 1.2.8 (2023-10-16)

### Bug Fixes

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Removed series name from fitSeriesProperties.
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed proportional confidence interval rendering.

## 1.2.7 (2023-09-11)

### Features

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Property panel changes:
  * Added tags on the dataframe level
  * Added proper tooltips on properties
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added errorModel property.

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added check for logarithm in the changeBounds method
  * Fixed confidence intervals for logarithmic modes

## 1.2.6 (2023-09-01)

### Features

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Stored the chart options in the dataframe tags
* T-Test now throws an error if sample size is less than or equal to 1.

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added checks for axes values

## 1.2.5 (2023-08-25)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added lineStyle to the IFitSeriesOptions interface
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Added code comments about line style
  
### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed the curve fitting algorithm
  * Fixed axes scales

## 1.2.4 (2023-08-23)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added standard deviation to the IFitPoint interface
  * Added marker size to the IFitPoint interface

## 1.2.3 (2023-08-12)

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
