# statistics changelog

## 1.12.11 (2026-08-28)

* MPO: Rearranged the profile editor internals — the operations behind its controls and the desirability curve math are shared API now

## 1.12.10 (2026-08-12)

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Fit: Added `IFitChartOptions.labels` for values that describe the whole plot and `showLabels` for the names to render. `IFitSeriesOptions.labels` keeps carrying per-curve values. Removed `labelOptions` and `IFitChartLabelOptions`, which nothing ever read - its `visible` is now `showLabels` and a label takes the colour of the series it belongs to. Stored JSON carrying the key is unaffected, since it never had an effect
* GROK-17637: Fit: Added `IFitChartData.explicit`, naming the options a user set at that level. Only an explicitly set level outranks the value a series declares for itself; anything absent stays advisory, as before
* GROK-17637: Fit: `FitConstants.TAG_FIT` moved from `.fit` to `.%fit`, the prefix the platform serializes into a layout, so column and dataframe level options survive a reopened project. The previous name is kept as `TAG_FIT_LEGACY`, read once and then removed, so only one of the two ever holds the options
* GROK-17637: Fit: Exported the typed fit results (`SigmoidFit`, `LinearFit`, `FourPLDoseResponseFit`, etc.), which expose named parameters instead of positional indices, plus `interceptY` and the raw `parameters`
* GROK-17637: Fit: Added `FitFunction.statisticsProperties`, describing the statistics each fit function produces
* GROK-17637: Fit: Added precomputed data points and log options to `fillParams`, and widened its series argument to `IFitSeries`
* GROK-17637: Fit: Fixed custom JS fit functions losing their fitted parameters, and `4pl-dose-response` reporting `ec50` instead of `ic50`
* GROK-17637: Fit: A custom JS fit function reports `interceptY` again, derived from the fitted curve as the pre-typed API did
* GROK-17637: Fit: Added `getSeriesFit`, `getStatistic` and `getStatisticProperty` - statistics are now derived from the typed fit, with legacy `FitStatistics` names resolved through an alias table
* GROK-17637: Fit: Added `toDataSpace`, converting log-fitted statistics back to data space in one place for both axes; `pIC50` is derived there, assuming a molar concentration
* GROK-17637: Fit: Added typed `getFitFunction` and `isFit`, and derived the statistic-name unions from the fit classes
* GROK-17637: Fit: Statistic labels now come from each fit function's `parameterNames`, and `maxY`/`minY` are derived from the asymptotes instead of mislabelling parameter slots
* GROK-17637: Fit: Tightened `IFitSeriesOptions` - removed the catch-all index signature and typed `fitFunction`, `showPoints` and `droplines`
* GROK-17637: Fit: Fixed `FitSeries.fit` being a circular own property that made `JSON.stringify` of a series throw
* GROK-17637: Fit: Broke the `fit-data` / `fit-engine` circular import by adding the `fit-points` leaf module and moving `fitSeries` next to the optimizer
* GROK-17637: Fit: Moved `LogOptions` to `fit-curve` alongside the other shared types. `getDataPoints`, `getMedian`, `getMedianPoints`, `logIC50ParameterBounds` and `fitSeries` moved out of `fit-data` into the modules that own them, so deep imports of those names from `fit-data` have to be updated
* GROK-17637: Fit: Removed the unused `getInvertedFunctions`, `FitInvertedFunctions` and `getInvError`, which inverted the sigmoid only via positional parameter indices
* GROK-17637: Fit: Removed the unused `FitFunctions` lookup class and the `FitSeries.fit` getter, both superseded by `getFitFunction` and `getSeriesFit`
* GROK-17637: Fit: `FitFunction` now builds `statisticsProperties` once in the base class from `statisticFields`, instead of each fit function repeating the memo
* GROK-17637: Fit: `statisticsProperties` no longer labels `top`/`bottom` as "Max Y"/"Min Y", which now name the derived asymptote statistics
* GROK-17637: Fit: Renamed `new-fit-API.ts` to `fit-engine.ts`, so the module is named for what it holds rather than for when it was added
* GROK-17637: Fit: A legacy statistic name a fit function does not produce now falls back to the parameter slot the pre-typed API read, so a `top` column extracted from a linear or exponential curve keeps its value instead of turning null. `getStatisticProperty` falls back the same way, so such a statistic still renders on the plot
* GROK-17637: Fit: `getSeriesFit` and `getSeriesConfidenceInterval` no longer write the fitted parameters back onto the series - they are in fit space, and the series contract is data space, so the next caller converted them twice
* MPO: Re-exported `desirabilityScore` for callers that score a single value (it was inlined into the column-at-a-time path and dropped from the public surface, breaking the PowerGrid build)
* MPO: Add ability to rename the property in the design mode
* MPO: A numerical property's min/max now default to its bound column's actual range instead of 0-1, locking once the user sets it by hand or saves the profile

## 1.12.9 (2026-07-16)

* Compute functions dialog: Optional per-function "re-run on open" checkbox (`rerunOnOpenOption`), with the flag surfaced in the dialog result and template compute types.
* Compute functions: `joinQueryResults` now overwrites an existing result column in place (matched by molecule) instead of appending a duplicate `<name> (2)` column.

## 1.12.8 (2026-07-03)

* MPO: Format score columns adaptively
* MPO: Add an option to invert numerical desirability profiles
* MPO: Add logarithmic scale support for numerical desirability profiles
* MPO: Fix a freeze when reopening dialogs containing large profiles

## 1.12.7 (2026-06-04)

* MPO: Optimized score calculation for large datasets

## 1.12.6 (2026-04-30)

* GROK-20056: MPO: Prevent rerendering on column input hover and wheel scroll

## 1.12.5 (2026-04-15)

* MPO: Use column input for auto-updating column selectors
* MPO: Fix column selector hover rerendering and layout overflow
* Compute dialog: Add single-select radio button mode for MPO functionality

## 1.12.4 (2026-04-08)

* MPO: Share compute function engine with HT

## 1.12.3 (2026-04-02)

* Statistics: Update fit consts

## 1.12.2 (2026-03-26)

* MPO: Add profile versioning

## 1.12.1 (2026-03-19)

* MPO: Added resize support for the line editor
* EDA: Fixed MPO line editor lookup issue in newer statistics library versions

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
