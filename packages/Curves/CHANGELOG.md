# Curves changelog

## v.next

* Flow: Added Add Curve Statistic and Add Aggregated Curve Statistic, which take a real `column {semType: fit}` slot. The `addStatisticsColumn` / `addAggrStatisticsColumn` originals address the curve column by NAME — a free-text field with no picker and no filter — and stay as they are because the Fit pane and the Data to Curves pipeline call them by name as replayable transforms
* Curve statistics: The free-text `propName` and `aggrType` parameters now declare `choices:` and a default, on both the originals and the new twins
* Fit Dose-Response Curves: Constrained the column pickers by type (numerical for concentration/readout, categorical for the identifiers) and marked the seven mandatory columns `nullable: false` — a column parameter defaults to nullable, which read as optional and let a half-configured node run
* Calculate MSR: Added captions, descriptions, and column-type constraints to the five unconstrained mandatory column inputs
* GROK-20232: Fit renderer: degrade gracefully on malformed fit-cell JSON instead of throwing on every grid redraw
* Fit options: Fixed a Column or Dataframe level change having no effect on a cell where the same property had been set at cell level. Whether any cell overrides a property was detected once and never revisited, so a cell that gained an override afterwards was invisible to the pass that clears them. The detection is now stamped with the column version it was computed from, and every cell write bumps that version, so it re-derives itself
* Fit renderer: Fixed the legend drawing each label over its own marker. The legend positions text by measured width, but did not set its own `textAlign`, so it inherited the centred alignment the axis labels leave behind and every label shifted half its width left
* GROK-17637: Statistics: Changing a chart or series option at Column or Dataframe level no longer recalculates the extracted statistic columns unless the option can affect a statistic - `mergeSeries` does not, since the renderer merges on a copy and the statistics read the original series. Clearing the per-cell overrides is how the option takes effect, but it was itself announced as a data change, so a colour or a title refit every dependent column
* GROK-17637: Statistics: Fixed an aggregated `interceptX` being reported in log space while `ic50` was a concentration - it is the default statistic of both aggregation transforms, so a saved project replayed -6 into a column of 1e-6 values, and an outlier toggle mixed both scales into one column
* GROK-17637: Statistics: The y asymptotes are reported in the space they were fitted in, as before. A stored `bottom` of 0 has no finite logarithm, so a forward and inverse conversion could never be exact inverses - an unguarded inverse reported that bottom as 1
* GROK-17637: Statistics: A merged-series cell no longer caches its fit under the key the first real series uses, which made an extracted statistic depend on whether the row had been painted
* GROK-17637: Statistics: `Add Aggregated Curve Statistic` no longer offers aggregations that are not converted back to data space, matching the other two aggregation functions
* GROK-17637: Statistics: A legacy statistic a fit function has no field for (`top` on a linear fit) now reads the same aggregated as it does per series, instead of coming back empty
* GROK-17637: Tests: Guarded the follow-on fixes - the 3DX parser not inventing `logY`, the detached recalculation path stripping the stray `|`, `interceptY` resolving to a descriptor on a custom fit, and the transform functions offering only aggregations that convert back to data space
* GROK-17637: Statistics: The renderer no longer replaces the series of the cached chart data when merging series, which made an extracted statistic depend on whether the row had been painted
* GROK-17637: Statistics: The recalculation path strips the stray `|` from a cell value like the initial parse does, so a row no longer blanks only when it is edited
* GROK-17637: Statistics: `interceptY` has a descriptor for every fit function, so a custom JS fit renders it instead of holding a value the plot skips
* GROK-17637: XML 3DX: `logY` is emitted only when the export declares it - always emitting the key stopped a Column or Dataframe level option cascading onto the cell
* GROK-17637: Statistics: Aggregation choices are limited to the ones that convert back to data space, on the transform functions as well as in the panel, and the new `curveStatistic`/`curveAggrStatistic` declare `semType: fit` and non-nullable inputs like their legacy twins
* GROK-17637: Tests: Added coverage for the renderer and the property panel - a repaint leaving the cached parameters in data space, values-changed driving statistic-column recalculation, the panel constructing for a saved legacy statistic, and legacy statistic names mapping onto the current ones
* GROK-17637: XML 3DX: Fixed boolean attributes being read as strings - `drawLine="False"` and `logX="false"` both came back true, so a converted curve drew a fit line the source asks to hide, and would have been fitted in log space on a linear axis
* GROK-17637: XML 3DX: A declared fit function that is not registered now falls back to sigmoid instead of being passed through as an unknown name
* GROK-17637: XML 3DX: `logY` is read from the chart settings, where it was previously ignored
* GROK-17637: Statistics: Plot and property panel now show the statistics the series' own fit function produces, instead of a fixed list that rendered NaN for non-sigmoid fits
* GROK-17637: Behaviour change: Statistics are labelled with each fit function's own parameter names - `AUC`, `R²`, and `Max`/`Hill`/`IC50`/`Min` for dose-response. `Max Y` and `Min Y` now name the derived asymptote statistics rather than `top`/`bottom`
* GROK-17637: Behaviour change: `slope` on a linear or log-linear fit now reports the slope; it previously reported the intercept. `top` on those fits keeps its old value, so `top` and `slope` report the same number
* GROK-17637: Statistics: IC50 is reported as a concentration everywhere, so the plot and the property panel no longer disagree
* GROK-17637: Behaviour change: IC50 is converted out of log space whenever the x axis is logarithmic. The conversion was previously skipped for stored values >= 1, and for `4pl-dose-response` entirely, so nM/uM-scale curves and every dose-response curve now report a different number - by orders of magnitude in saved projects
* GROK-17637: Behaviour change: Statistics aggregated across series are averaged in fit space, so an averaged IC50 is now a geometric mean - 1e-7 and 1e-5 give 1e-6 rather than 5.05e-6. Counts, sums and shape statistics stay on their own scale
* GROK-17637: Added pIC50, derived as -log10(IC50). It assumes the concentration is molar, so it reads 6 too low on uM data and 9 too low on nM data
* GROK-17637: Options: A custom JS fit function is now named in the property panel instead of rendering a blank input
* GROK-17637: Options: Fixed the statistics and fit-function inputs binding to a throwaway copy of the options, so changing either did nothing
* GROK-17637: Statistics: Toggling an outlier refreshes extracted statistic columns by name off the typed fit, so columns named beyond the seven legacy statistics (`ic50`, `pIC50`, `maxY`) stay populated
* GROK-17637: Statistics: Fixed a recalculated statistic ignoring column-level chart options - a column-level `logX`, `allowXZeroes` or fit function was lost on the detached recalculation column
* GROK-17637: Statistics: A cell mixing fit functions now aggregates only the statistics every series produces, instead of averaging one over the subset that has it
* GROK-17637: Statistics: Changing a chart option at Column or Dataframe level now refreshes extracted statistic columns - the option lives in a tag, so nothing marked the curve column changed and the columns kept stale numbers while the plot updated
* GROK-17637: Statistics: Statistics that only carry meaning in data space (IC50/EC50 under a logarithmic x axis, and pIC50) are no longer offered for aggregations that are not converted back, where they would have been reported in log space or come back empty
* GROK-17637: Fixed the series colour type being widened to `string`, which disabled colour-option checking
* GROK-17637: Statistics: Added `curveStatistic` and `curveAggrStatistic` - the property panel now adds a statistic as a calculated column that recalculates when its curve changes and is recorded in the table's creation script
* GROK-17637: Statistics: Aggregation now covers every statistic the cell's fit functions produce rather than a fixed list, and the aggregated pane uses the same per-fit-function descriptors as the single-series one
* GROK-17637: Fixed confidence intervals reading a non-existent `series.dataPoints`, which discarded the cached points and recomputed them on every frame
* GROK-17637: Statistics: Moved the statistics calculation out of the property-panel handler and the package entry point into `fit/fit-statistics.ts`
* GROK-17637: Extracted the chart-data assembly and fit caches into `fit/fit-chart-data.ts`, breaking the `fit-renderer` <-> `fit-statistics` import cycle
* GROK-17637: Statistics: Fixed the reported IC50 of nM/uM-scale curves drifting (100 -> 2 -> 0.301). `IFitSeries.parameters` are now always in data space: neither the renderer nor the statistics convert them in place or write fitted parameters back onto the cached series

## 1.12.0 (2026-04-02)

### Features

* Added PZFX (GraphPad Prism) file viewer and handler
* Added dynamically extendable format support

## 1.11.0 (2026-03-13)

### Features

* Added editable GridCellWidget support

## 1.10.15 (2025-12-02)

* Png Renderer: Fix flickering and moved to PowerGrid

## 1.10.12 (2025-09-16)

* GROK-18808: Provided consistent names for Data top menu
* Plates: Skip rendering too small plates

## 1.10.11 (2025-09-12)

* Data-To-Curves: Well level additional columns

## 1.10.10 (2025-09-09)

* Extended PNG detector

## 1.10.9 (2025-08-27)

* MultiCurveViewer: Added showOutliers parameter

## 1.10.8 (2025-08-21)

* Data-to-Curves: Correct exclude heuristic

## 1.10.7 (2025-08-21)

* Rendering: Added showOutliers parameter

## 1.10.6 (2025-08-20)

* TestData: Wells: Improved concentration/activity, samples, barcodes generation for plates
* TreeBrowser function: Removed unnecessary browsePanel input
* Added MSR script
* Updated MSR env

## 1.10.5 (2025-08-20)

## 1.10.4 (2025-08-14)

* MultiCurveViewer:
  * Fixed issues with min/max x/y
  * Fixed slider is added on min/max x/y fields
  * Added legendColumnName property

## 1.10.3 (2025-08-01)

## 1.10.2 (2025-07-31)

* Data-to-curves: Added joining options

## 1.10.1 (2025-07-30)

* Data-to-curves: Fixed mismatched keys

## 1.10.0 (2025-07-29)

### Features
* 4plDoseResponse: Added 4PL dose-response curve fitting function
* Data-to-curves: fully covered script with history, datasync and 2 tier support
* Raw PNG renderer

## 1.8.0 (2025-03-31)

### Features

* GROK-17295: Fitting caching:
   * Performance improvements
   * Cached fitting and wrote tests
   * Fixed margins
* Fit: Disabled caching for non valid row indexes
* Outliers: Made outlier toggling not to invalidate the whole grid
* Added custom event for outlier toggling
* Added 4pl-regression fit function
* #3195: Plate Viewer:
   * Supported plate handling
   * Implemented automatic plate import from Excel
   * Added new fluent API
   * Enabled DRC with grouping
   * Supported multi-plate file folders
   * Added details, stats column, IC50, curve preview, statistics recalculation, outlier marking, configuration options
   * Implemented content-dependent file import
   * Added rendering for 1536-well plates
   * Enabled inspection of single series in plate
   * Added demo
* Plate readers: Added Delfia Envision and Spectramax parsers
* Added docs on plate DRC

### Bug fixes

* Fixed MultiCurveViewer issues
* Fit panel: Fixed namings and inputs
* GROK-17492: Fixed error when adding statistics column if the series has no name

## 1.6.0 (2024-09-09)

* Fixed tests

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
