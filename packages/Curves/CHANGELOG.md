# Curves changelog

## 1.5.0 (2024-08-23)

### Features

* Improved fitting algorithm - didn't fit properly because was stuck in the bounds

## 1.4.6 (2024-08-13)

### Features

* Added exponential fit function

## 1.4.5 (2024-07-24)

## 1.4.4 (2024-07-02)

### Bug fixes

* [#2924](https://github.com/datagrok-ai/public/issues/2924): Fixed connectDots doesn't work if fit is turned off

## 1.4.3 (2024-05-23)

### Bug fixes

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Changed labels (Marker Type -> Marker and Outlier Marker Type -> Outlier Marker)

## 1.4.2 (2024-05-21)

### Features

* [#2855](https://github.com/datagrok-ai/public/issues/2855): Added outlier marker option

## 1.4.1 (2024-05-21)

### Features

* [#2797](https://github.com/datagrok-ai/public/issues/2797): Added log-linear fit function

## 1.4.0 (2024-05-11)

### Features

* Changed optimizer for fitting

### Bug Fixes

* [#1645](https://github.com/datagrok-ai/public/issues/1645): MultiCurveViewer: Removed unnecessary margins

## 1.3.1 (2024-05-03)

### Features

* [#1645](https://github.com/datagrok-ai/public/issues/1645): MultiCurveViewer:
  * Added filter on fit columns in the MultiCurveViewer
  * Made property panel fields work and made automatic chart options merge
  * Code refactoring and restructure, render improvements

## 1.3.0 (2024-04-18)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improve curves properties and rendering:
  * Added tooltips on points
  * Added icon to MultiCurveViewer
  * Added mergeSeries property and method
  * Added legend rendering with column labels
* [#2754](https://github.com/datagrok-ai/public/issues/2754): Implemented capability just to connect the points (without fitting)
* [#1645](https://github.com/datagrok-ai/public/issues/1645): MultiCurveViewer:
  * Added columns selection
  * Added column series merging
  * Added colors and fixed error messages
  * Added curves limit and fixed showColumnLabel
  * Turned off droplines and confidence intervals
  * Made same fit line style in case of 20+ curves

### Bug Fixes

* [#2748](https://github.com/datagrok-ai/public/issues/2748):
  * Fixed whole table is broken if cell with curves contains specific data
  * Added red cross if the curve is broken
  * Added error text

## 1.2.16 (2024-03-15)

### Features

* [#2105](https://github.com/datagrok-ai/public/issues/2105): Added outlier color property.

### Bug Fixes

* Renamed "Show curve confidence intervals" to "Confidence intervals"
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed statistics calculation
* [#2103](https://github.com/datagrok-ai/public/issues/2103): Made statistics columns adding functions DG functions

## 1.2.15 (2023-11-17)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added capability to substitute zeroes in curve fitting for logarithmic data

### Bug Fixes

* Fixed MultiCurveViewer throws an error when adding to an arbitrary dataset

## 1.2.14 (2023-11-10)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Changed droplines min size for rendering

## 1.2.13 (2023-11-09)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed curves not rendering because of undefined x if branch.
  * Fixed linear fit not rendering.

## 1.2.12 (2023-11-02)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed exception if no x or y coordinates present.
* [#2103](https://github.com/datagrok-ai/public/issues/2103): Property panel changes:
  * Returned null values if lack of data in aggregated columns
  * Replaced DG.Stats.fromColumn() with DG.Stats.fromValues()

## 1.2.11 (2023-10-26)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed exception thrown on adding Form viewer.

## 1.2.10 (2023-10-16)

### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed proportional confidence interval rendering.

## 1.2.9 (2023-10-11)

### Features

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Implemented aggregations for series statistics
  
### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Fixed curves if no points present in series
* [#2103](https://github.com/datagrok-ai/public/issues/2103): Property panel changes:
  * Fixed the color input for the curves
  * Inserted stat column next to curves column
  * Fixed parameter column adding
  * Fixed adding statistics column for null undefined series
  * Replaced the MultiCurveViewer with CellRenderViewer in the property panel
  * Fixed statistics calculation and rendering
* [#2394](https://github.com/datagrok-ai/public/issues/2394): Replaced autostart tag with init tag

## 1.2.8 (2023-09-25)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved the curves demo app
  
### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Don't allow clickToToggle if small cell size
  * Fixed margins in small cells

## 1.2.7 (2023-09-11)

### Features

* [#2103](https://github.com/datagrok-ai/public/issues/2103): Property panel changes:
  * Added color for the series in the property panel
  * Added capability to store chart options on the dataframe level in tags
  * Restructured the property panel
  * Added dataframe and plot-only switch modes in the property panel
  * Made the accordion to restore the property panel state
  * Added proper tooltips on properties
  * Implemented the capability to switch the plot state on the property panel
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Decreased the sizes for the plot title rendering
  * Changed the min axes cell width
  * Added errorModel property
* [#2100](https://github.com/datagrok-ai/public/issues/2100): Made cell renderer edit mode a resizeable window
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Added documentation about property panel and error model
  
### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added dataBounds checks in case of logarithmic values
  * Fixed droplines rendering for linear and logarithmic modes
  * Fixed the tooltip for small cells
  * Fixed confidence intervals rendering for logarithmic modes

## 1.2.6 (2023-08-25)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Added different line styles rendering
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Added documnetation about line styles
  
### Bug Fixes

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Fixed the curve fitting algorithm
  * Fixed axes scales

## 1.2.5 (2023-08-23)

### Features

* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering:
  * Added standard deviation rendering for points
  * Added custom point color coding
  * Added custom point marker
  * Added custom marker size rendering
* [#2106](https://github.com/datagrok-ai/public/issues/2106): Updated the documentation regarding the standard deviation, the point color, marker type and size

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
