# Charts changelog

## 1.5.0 (WIP)

### Features

* Added Multiplot viewer

## 1.4.4 (2024-02-27)

### Features

* [#3249](https://github.com/datagrok-ai/public/issues/3249): Tree: Add cross-viewer selection

### Bug Fixes

* [3245](https://github.com/datagrok-ai/public/issues/3245): Tree: Usability issues:
  * Combination "Row Source: All - On Click: Filter" works differently in several cases
  * The viewer re-renders on every branch click, even when the row source remains unchanged
  * Pressing ESC on the keyboard deselects rows in the table, but the highlighting remains on the viewer
* GROK-17619: Charts | Tree viewer: Symbol size changes when aggregation function is set

## 1.4.3 (2024-01-30)

### Features

* [#3221](https://github.com/datagrok-ai/public/issues/3221): Tree: Improvements:
  * Display row counts for each branch
  * Allow selection of sets in each branch to be either filtered or selected in the grid
  * Enable font size customization
  * Fix issue where the last-level categories are not readable in horizontal view
  * Add a setting to hide null values
  * Support structure rendering
  * Use standard tooltips

### Bug Fixes

* GROK-17376: Charts | Tree: "NaN.floor()" error while changing Color in properties
* GROK-17405: Charts | Tree: Changing Size in properties caused errors
* GROK-17414: Charts | Tree: Dragging down the three causes errors
* GROK-17420: Charts | Tree: Prevent view reset when changing Symbol size, Label rotation, or Font size
* GROK-17424: Charts | Tree: Flickering structures when ShowMouseOverLine is enabled

## 1.4.2 (2024-12-19)

### Bug Fixes

* GROK-16954: Charts: Fixed globe throws errors if no column found
* Charts: Sunburst: Improved logic for labels rendering
* Charts: Word cloud: Fixed cases when no column selected
* [#3196](https://github.com/datagrok-ai/public/issues/3196): Charts: Sunburst causes performance issues
* Charts: Sunburst: Resolved an issue where certain columns were not used when applying the row source
* Charts: Sunburst: Increased the size of rendered structures in the tooltip and implemented dynamic size adjustments based on segment size

## 1.4.1 (2024-10-03)

### Bug Fixes

* GROK-16844: Word cloud: Fixed viewer opens empty on SPGI dataset
* Charts: Sunburst: Fixed test throws exception when this.dataframe is null in _render() in then
* Charts: Tree: Fixed test throws exception when this.dataframe is null in _render() in then

### Features

* [#3090](https://github.com/datagrok-ai/public/issues/3090): Sunburst: Usability improvements:
  * Reduce free space around
  * Reordering of the columns in the dialog is not reflected on the fly
* [#3097](https://github.com/datagrok-ai/public/issues/3097): Sunburst: change Row Source to 'All' once the setting 'On Click' is switched to 'Filtered'
* Charts: Sunburst: Improve logic for labels rendering

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
