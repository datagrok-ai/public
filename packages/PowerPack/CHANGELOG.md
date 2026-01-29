# Power Pack changelog

## 1.7.19 (2025-12-19)

* DBExplorer: Handle semtypes

## 1.7.15 (2025-10-27)

* Removed recent projects widget

## 1.7.14 (2025-10-20)

* GROK-19080 Add new column: type change is broken

## 1.7.12 (2025-10-15)

* Add new column dialog: Some validation fixes

## 1.7.11 (2025-10-15)

* Add new column dialog: Case insensitive packages search for hints

## 1.7.10 (2025-10-09)

* Add new column dialog: Fixed validation error when passing empty array to vector functions

## 1.7.9 (2025-10-08)

* GROK-18237: Widgets: Removed the Settings button for widgets
* GROK-18139: Activity dashboard:
  * Made sure clicking on entities will open the property panel (needs new JS API version)
  * Removed `Previous tip` arrow
  * Retired `Learn` widget and added it to the activity dashboard instead

## 1.7.8 (2025-10-07)

* GROK-18139: Activity dashboard: Made subwidget functions return DG.Widget instead of HTMLElement[]
* Removed unnecessary logs
* Removed DG User handler

## 1.7.7 (2025-10-02)

* GROK-18139: Activity dashboard:
  * Fixed tips of the day
  * Removed admin tab
  * Added check for new users
* GROK-17174: Customize widgets: Fixed menu is empty
* GROk-18533: Browse: Fixed scroll position is reset after opening demo

## 1.7.6 (2025-10-02)

* GROK-18139: Activity dashboard: Fixed getting recent entities doesn't work correctly
* GROK-18977: Activity dashboard: Fixed text overlapping on window reduced width

## 1.7.5 (2025-09-26)

* GROK-18139: Activity dashboard: Fixed styles, fixed getting demo app hierarchy

## 1.7.4 (2025-09-25)

* GROK-18139: Activity dashboard:
  * Added ui.wait
  * Optimized performance
  * Merged requests to one where applicable
  * Fixed header jumps when removed
  * Refactored code
  * Improved styles
  * Added scroll in case something left below
  * Removed unnecessary data connections
  * Fixed notifications
  * Added link to reports
  * Added content for new users
  * Made sure content from other package funcs also comes to the dashboard
  * Uncommented activity dashboard
  * Moved tip of the day for all users
  * Fixed errors
  * Added null check
  * Fixed tips of the day
  * Fixed styles
  * Added randomized demo of the day
  * Added Demo Datasets subwidget, fixed it for new users
  * Fixed subwidgets names
  * Moved random tip to the left and to the bottom
  * Removed unnecessary demo datasets, added demos/tutorials
  * Fixed styles for the tip of the day
  * Added capability to change tips
  * Added prev/next steps and tutorials of the day
  * Added generalized tips of the day based on day and week (same for all the stands)
* Scatter plot: Adjacent column not changing upon changing the axes in formula line preview
* Bump db-exp
* GROK-18733: Recent projects:
  * Added ui.wait()
  * Made sure ui.wait() located in the center
* Top menu: Removed formula lines
* Excel: Improved cases when same column headers in the file

## 1.7.3 (2025-08-22)

* Commented out activity dashboard

## 1.7.2 (2025-08-11)

* GROK-18656: Power Search: Error when pressing ENTER in some cases
* Excel:
  * Fixed tests timeouts
  * Fixed sometimes empty values break the parsing
* GROK-18139: Activity dashboard:
  * Added scroll on each subwidget if needed
  * Removed the header

## 1.7.1 (2025-08-11)

## 1.7.0 (2025-08-07)

* DBEXP: Data source name config
* Return string error instead of throwing error in validation function
* md/mdx file preview fixes
* GROK-17445: Complex calculated columns persisting with the layout
* Added ability to disable viewers in the viewers gallery
* GROK-17282: Fixed cases when formula contains no column after = but 2 columns before it
* Add new column:
  * Added validation for typed lists like num_list, string_list etc.
  * Skip dataframe param validation for vector functions
  * Small UI fixes
  * Additional column type validation, fixed linter errors
  * Warning in case renaming to existing name
* [#2449](https://github.com/datagrok-ai/public/issues/2449): Add new column: Using non-empty rows in preview
* Renamed disableInViewerGallery tag to showInGallery
* [#3306](https://github.com/datagrok-ai/public/issues/3306): Calculated columns: If function types validation error
* File handlers: Added ExcelJS xlsx file handler with web workers, added tests
* [#3389](https://github.com/datagrok-ai/public/issues/3389): Formula lines: Fixed formula lines preview for Scatterplot used in Trellis doesn't reflect axes configured in Scatterplot
* [#3269](https://github.com/datagrok-ai/public/issues/3269): Projects: Added sheetName for JS file-handler (Excel in this case), removed unnecessary code
* Edit formula widget improvements
* [#3392](https://github.com/datagrok-ai/public/issues/3392): Scatterplot: Fixed axes in formula lines dialog should be the same as on the viewer
* GROK-18355: Add new column: Fixed receiving errors when trying to type name
* Renamed formula widget
* GROK-18377: Line Chart: Fixed opening the Formula Lines dialog causes an error if the X-axis is split by time
* GROK-18397: Line Chart: Fixed empty formula lines dialog for some splits enabled when date/time columns are set for X/Y
* GROK-18398 : Complex calculated columns: unexpected errors in Editing
* GROK-18396: Line Chart: Formula lines are not displayed in the dialog when the X-axis is time-splitted fixed
* Search improvement
* App and samples search
* Lazy loaded object handler
* Removed denial search
* GROK-18426: Add new column: Fixed renaming column to an existing name throws errors
* GROK-14990: Viewers gallery: Made core viewers disabled if they can't be added in the gallery
* GROK-18426: Add new column: Fixed renaming column to an existing name throws errors
* GROK-18139: Activity dashboard: Implemented and added Activity dashboard

## 1.5.0 (2025-03-29)

## 1.4.1 (2025-03-18)

### Bug Fixes

* GROK-17755: Failed to open "Add viewer" dialog

## 1.4.0 (2025-02-24)

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
