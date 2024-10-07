# Charts changelog

## 1.4.1 (2024-10-03)

### Features

* [#3090](https://github.com/datagrok-ai/public/issues/3090): Sunburst: Usability improvements:
  * Reduce free space around
  * Reordering of the columns in the dialog is not reflected on the fly

### Bug Fixes

* Dependency on the chem package has been removed

## 1.3.5 (2024-08-30)

### Features

* Charts: Sunburst: Speed up

## 1.3.4 (2024-08-23)

### Bug Fixes

* [#2994](https://github.com/datagrok-ai/public/issues/2994): Charts | Sunburst: Double click on the viewer resets the view
* [#2992](https://github.com/datagrok-ai/public/issues/2992): Charts: Sunburst does not select/filter on empty category
* Charts: Sunburst: Update emphasis in order not to fade out other data
* Charts: Sunburst: Use only categorical columns for building the hierarchy
* [#2979](https://github.com/datagrok-ai/public/issues/2979): Charts: Sunburst is showing default columns instead of selected ones after applying layaout in some cases
* Charts: Demo: Fix globe viewer doesn't show points

## 1.3.3 (2024-08-08)

### Bug Fixes

* [#2954](https://github.com/datagrok-ai/public/issues/2954): Sunburst: Unsupported columns are not filtered out in column selector

## 1.3.2 (2024-07-24)

## 1.3.1 (2024-07-02)

### Features

* GROK-15340: Charts | Radar: Support multiple series

### Bug Fixes

* GROK-15132: Charts | Timelines: Error when using the legend

## 1.3.0 (2024-03-21)

### Features

* [#2500](https://github.com/datagrok-ai/public/issues/2500): Sunburst: Improvements:
  * Implemented on-click filtering
  * Added Row Source
  * Provide a dialog like Order and Hide with all/visible/ hidden
  * Turned off the ability to drill down on click (ctrl +click or cmd+click)
  * Implemented the logic of wrapping the long labels
  * Added a tooltip: label + row count

### Bug Fixes

* Radar viewer minor imrpovements
* Links fixes
* [#2500](https://github.com/datagrok-ai/public/issues/2500): Sunburst: Labels were hidden after the structures were added
* GROK-15166: Charts | Sunburst: Error when trying to switch the Table property

## 1.2.2 (2023-09-19)

### Bug Fixes

* Fixed error when closing Radar viewer
* Fixed errors in Surface plot
* [#2096](https://github.com/datagrok-ai/public/issues/2096): Sunburst: Added different click types (make work for Mac keyboard shortcuts)

## 1.2.1 (2023-08-24)

### Bug Fixes

* [#2292](https://github.com/datagrok-ai/public/issues/2292): Suburst: Doesn't respond on table switching

## 1.2.0 (2023-08-17)

### Bug Fixes

* [#2098](https://github.com/datagrok-ai/public/issues/2098): Sunburst: UI/UX improvements:
  * Fixed an error when clicking on the root node
  * Minor tooltip and rendering fixes
  * Render if column semtype is molecule

## 1.1.1 (2023-08-02)

### Features

* [#2098](https://github.com/datagrok-ai/public/issues/2098): Sunburst: Added group tooltips using showRowGroup

## 1.1.0 (2023-08-01)

This release focuses on adding new functionality and improving the existing one.

### Features

* Switched to class decorators for registering viewers
* [#1974](https://github.com/datagrok-ai/public/issues/1974): Timelines: Set the default date format based on the time difference between event start and end values
* [#2096](https://github.com/datagrok-ai/public/issues/2096): Sunburst: Add different click types, highlighting, selection, remove animations:
  * Implemented Ctrl + Click for selecting specific parts of the plot
  * Implemented Shift + Click for adding selection on the plot
  * Implemented Ctrl + Shift + Click for removing selection
  * Implemented highlighting to sync with other viewers
  * Removed animations
* [#2098](https://github.com/datagrok-ai/public/issues/2098): Sunburst: UI/UX improvements:
  * Added tooltips
  * Added molecule rendering on plot
  * Synchronized plot with the table color coding
  * Text doesn't render if it doesn't fit
* [#2097](https://github.com/datagrok-ai/public/issues/2097): Sunburst: Add reset viewer functionality:
  * Added viewer reset on double click and to the context menu
  * Unsubscribes from events when viewer detached

### Bug Fixes

* [#1968](https://github.com/datagrok-ai/public/issues/1968): RadarViewer does not allow changing columns in SpiderChart:
  * Show actual number of columns that are selected
  * RadarViewer contains non numerical columns
* [#1962](https://github.com/datagrok-ai/public/issues/1962): RadarViewer errors on bigint datatype
* Surface plot fixes (format, axis labels)
* [#2098](https://github.com/datagrok-ai/public/issues/2098): Sunburst: UI/UX improvements:
  * Verified that the hierarchy column order is persisted properly if we reopen the "hierarchy" dialog
  * Fixed wrong behavior when second-level classes were not aggregated

## 1.0.24 (2023-05-26)

### Bug Fixes

* [#1973](https://github.com/datagrok-ai/public/issues/1973): Timelines: Visualization is not shown if the Split by column contains long category names
