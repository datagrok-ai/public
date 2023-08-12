# Curves changelog

## 1.2.4 (2023-08-12)

### Features

* [#2106](https://github.com/datagrok-ai/public/issues/2106): Updated the documentation about parameter order and plot title
* [#2105](https://github.com/datagrok-ai/public/issues/2105): Made Curves by default as a DG categorical color
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added x- and y-axes labels rendering
  * Added title rendering if the cell size is enough

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed log curves rendering
  * Fixed sizes for axes rendering

## 1.2.3 (2023-08-07)

### Bug Fixes

* [#2104](https://github.com/datagrok-ai/public/issues/2104): Fixed fit detector

## 1.2.2 (2023-08-04)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Changed minBound and maxBound to min and max in parameterBounds
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Updated the documentation about parameterBounds

## 1.2.1 (2023-08-02)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added dropline rendering for IC50
  * Added empty cell value handling
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Wrote TS docs and extended the documentation
* Added error handling in the JnJ parser

## 1.2.0 (2023-07-21)

This release focuses on adding new functionality and improving the existing one.

### Features

* Added outliers switching in grid cells.
* Added different modes of rendering: candlesticks (with box plot statistics) and both (with candlesticks outliers).
* Improved fit detector.
* Improved curves and confidence intervals rendering for better smooth.
* Added the capability to run curves by user-defined Javascript function with caching.
* Improved rendering in the small cells.
* Added different marker types.
* Improved axes cell rendering.
* Implemented interactive parameter recalculation.
* [#2105](https://github.com/datagrok-ai/public/issues/2105): Outliers are now shown in red color and with a much bigger size.
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Curves now don't render before the axes start.
  * Changed margins in grid cell for axes.
  * Added the capability to control parameter bounds.
* [#2100](https://github.com/datagrok-ai/public/issues/2100): Made cell renderer edit mode a resizeable window.

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed points rendering in logarithmic data.
  * Fixed curves rendering by inverting the y while converting the coordinates.

## 1.1.0 (2023-05-19)

This release focuses on improving the analysis stability and usability.

### Features

* Added support for log transform and axis inverse.
* Improved curve fitting and confidence intervals:
  * Added a semantic type detector for "fit" to render fitted curves.
  * Automated data sorting during fitting.

### Bug Fixes

* Fixed performance issues during fitting.
