# Power Pack changelog

## 1.4.1  (2025-03-18)

### Bug Fixes

* GROK-17755: Failed to open "Add viewer" dialog

## 1.4.0  (2025-02-24) - 1.3.0 (2024-12-24)

### Features

* [#2708](https://github.com/datagrok-ai/public/issues/2708): Formula lines: Changed "\$\{x\}" to "x" on Formula Lines
* Updated validation, added tests for BinBySpecificLimits
* Widgets visibility
* PowerPack: Add new column dialog: Added code mirror border when not focused
* Viewers gallery: Temporarily removed Shape map (until we reincarnate it)
* Scrolling: Saved scrolling in tab mode for browse and other things

### Bug Fixes

* GROK-17109: Calculated columns: Fixed columns are not saved to project with data sync
* Welcome: Widgets: not showing the "drag" cursor on the header since it's not working
* Fixed table opening from powerSearch
* Fixed widgetHost func
* Fixed error, when returned widget is null

## 1.2.3 (2024-11-13)

### Bug Fixes

* [#3146](https://github.com/datagrok-ai/public/issues/3146): Add new columns: Error on calculated columns validation if a complex formula is pasted
* [#3147](https://github.com/datagrok-ai/public/issues/3147): Add new columns: Columns are not calculated for some formulae inside if clause

## 1.2.2 (2024-10-28)

### Bug Fixes

* [#3119](https://github.com/datagrok-ai/public/issues/3119): Calculated column: formula with multi-argument functions is parsed incorrectly if more than one argument contains a column

## 1.2.1 (2024-09-27)

Datagrok DB explorer

## 1.2.0 (2024-09-09)

### Features

* Add new column autocomplete: show packages in the separate section in the bottom of list
* Add new column: Filter out irrelevant functions for autocomplete

### Bug Fixes

* GROK-16608: Add new column: Function 'if(true, "yes", error("Error"))' leads to validation error
* Add new column: Code mirror focus is removed after autocomplete
* [#3009](https://github.com/datagrok-ai/public/issues/3009): Calculated column hints break with round bracket + curly bracket
* Add new column: Fixed preview grid styles
* [#3008](https://github.com/datagrok-ai/public/issues/3008): Qnum function is missing from hints

## 1.1.13 (2024-08-27)

### Bug Fixes

* Fixed the PDB search (only showing existing proteins)
* Add new column: Some validation fixes

## 1.1.12 (2024-08-21)

### Bug Fixes

* [#2918](https://github.com/datagrok-ai/public/issues/2918): Typing $ in the calculated column dialog pops up the available column drop down also within a quoted field

## 1.1.11 (2024-08-20)

### Bug Fixes

* [#2993](https://github.com/datagrok-ai/public/issues/2993): Calculated columns dialog: Cannot type in any input after adding a calculated column
* [#2988](https://github.com/datagrok-ai/public/issues/2988): Calculated columns dialog: Fixed qnum type validation
* [#2981](https://github.com/datagrok-ai/public/issues/2981): Calculated columns dialog: Typing calculated formula: input is overwritten unexpectedly in some cases

## 1.1.10 (2024-08-07)

### Bug Fixes

* Add new column: incorrect input parameters when autocomplete

## 1.1.9 (2024-08-05)

### Bug Fixes

* [#2949](https://github.com/datagrok-ai/public/issues/2949): Calculated columns dialog: extra scroll bars

## 1.1.8 (2024-07-24)

* Dependency: datagarok-api >= 1.20.0*

### Features

* Updated 'Add new column' dialog: autocomplete and hints, formula validation, sorting by column type

## 1.1.7 (2024-06-12)

* Semantic value extractors
* Better tooltip phrasing
* [#2747](https://github.com/datagrok-ai/public/issues/2747): Formula lines:
  * Fixed vertical line is rendered as horizontal in the dialog
  * Added capability to change axes in line chart dialog
* Startup speed improvements
* Home page: Search: support for recognized entities
* [#2424](https://github.com/datagrok-ai/public/issues/2424): Viewers gallery:
  * Added invalid data handling in viewers gallery for Charts
  * Made Charts viewers disabled if the data is invalid in viewers gallery
* [#2830](https://github.com/datagrok-ai/public/issues/2830): Recent projects: Fixed double-clicking issue
* Help: Synchronized icons with the settings

## 1.1.5 (2023-11-28)

### Bug Fixes

* 1.18.0 Platform version compatibility, Recent projects widget

## 1.1.4 (2023-11-28)

### Bug Fixes

* Disabled queries search to increase performance, until DG 1.18.0

## 1.1.1 (2023-11-09)

### Bug Fixes

* [#2488](https://github.com/datagrok-ai/public/issues/2488): Fixed colour selector is missing for formula lines dialog

## 1.1.0 (2023-10-31)

### Features

* Provided welcome view
* Balloon is shown if viewer cannot be added
* Provided tooltip for long names
* Added toolbox toggle to statusbar

### Bug Fixes

* Fixed links
* Fixed wiki links
* Removed DistributionProfilerViewer

## 1.0.8 (2023-06-21)

Introduces Windows Manager, and a bunch of improvements in other areas

### Features

* Windows Manager (icons on the right of the status bar that control visibility of tool windows)
* Improved "Add New Column" dialog
* Redesigned Viewer Gallery
