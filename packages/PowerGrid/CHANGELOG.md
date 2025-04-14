# Power Grid changelog

## 1.5.5 (2025-03-07)

* GROK-17680: ImageUrl: Supported images with Datagrok paths, added docs

## 1.5.4 (2025-03-03) - 1.4.11

### Features

* Made pinned columns apply the same column header style as the standard grid
* GROK-17292: Implemented normalization for sparklines (row/column/global) with backward compatibility
* GROK-17359: Enabled renderer to define context values for sparklines in the context panel
* GROK-17600: Added zeroScale to sparklines, allowing them to be based at zero or a minimum value
* Allowed cross-origin images in canvas
* Tags cell renderer:
  * Added support for boolean columns
  * Minor rendering improvements
* Vlaaivis renderer:
  * Ensured no borders or filled segments for cells smaller than 40x40
  * Render only a light border

### Bug Fixes

* Fixed pinned columns not applying grid cell styles
* GROK-14428: Fixed duplicated row numbers in pinned columns
* GROK-16985: Fixed text not visible when editing in pinned columns
* [#3229](https://github.com/datagrok-ai/public/issues/3229): Fixed incorrect order in pinned columns after reordering by another column
* GROK-17599:
  * Different normalizations applied to all sparklines
  * Proper rendering in the context panel.
  * Fixed number formatting
  * Added margins and column name in the context panel
  * Fixed cases where `currentRow = -1`
* Fixed NaN errors in sparklines
* Fixed cases where the column is null in sparklines settings
* GROK-17550: Fixed error when saving to a project in some cases
* Fixed alignment in summary column tooltips
* GROK-17303: Fixed column name and value misalignment
* GROK-17293: Fixed invisible text in SmartForm when the column text is colored

## 1.4.10 (2024-11-18)

* Grid scrolling and current row fixes for pinned columns

## 1.4.5 (2024-10-25)

* WebGPU colors and different marker types for ScatterPlot

## 1.4.3 (2024-10-11)

* Fix webGPU Sc for non 4-divisible rows

## 1.4.2 (2024-10-11)

* Add webGPU ScatterPlot support

## 1.3.1 (2024-07-24)

### Features

* Introduced segments in PieChart

### Bug Fixes

* Added light gray circle border and fixed rendering when width is smaller than expected in PieChart
* Fixed name null in createTooltip in sparklines

## 1.3.0 (2024-02-20)

### Features

* [#2673](https://github.com/datagrok-ai/public/issues/2673): Implemented MultiChoice cell renderer.
* [#2681](https://github.com/datagrok-ai/public/issues/2681): Implemented Tags cell renderer.

## 1.2.6 (2023-11-30)

### Bug Fixes

* Forms viewer: Fixed indicator margins.

## 1.2.5 (2023-11-29)

### Features

* Forms viewer:
  * Updated viewer layout.
  * Added styles for updated molecular input.
  * Click on input makes row current.
  * Made column label clickable.
  * Added tooltip on value.
  * Added close button for column labels.
  * Added indicators for current and mouse over row.
  * Added default bottom position for Forms viewer.

### Bug Fixes

* [#2509](https://github.com/datagrok-ai/public/issues/2509): Fixed pinned columns in multiple views for the same table cause performance issues in some cases
* [#2528](https://github.com/datagrok-ai/public/issues/2528): Fixed row headers are lost after pinning and un-pinning columns in some cases
* Forms viewer: Fixed setting input text color.

## 1.2.4 (2023-11-27)

### Features

* Forms viewer:
  * Handled color coding changing, columns removing, input values changing.
  * Made inputs read-only.
  * Set default limit for columns number.
  * Added reaction to Ctrl, Shift, and Ctrl+Shift-click.
  * Added help link.
  * Updated viewer layout.
  * Added Forms viewer icon.

### Bug Fixes

* Forms viewer:
  * Fixed styles.
  * Fixed removing several columns at once.
  * Fixed freezing when mouse enters form.
  * Fixed rendering when showCurrentRow is unset.
  * Fixed rendering in case showCurrentRow and showMouseOverRow are both switched off.
  * Fixed text color coding.
  * Fixed showSelectedRows and checkboxes height.

## 1.2.3 (2023-11-24)

### Features

* Added image to context panel on current row.

### Bug Fixes

* Fixed double-click popup window.

## 1.2.2 (2023-11-20)

### Features

* Added Forms viewer.

### Bug Fixes

* Fixed molecules rendering in smartforms.

## 1.2.1 (2023-11-07)

### Bug Fixes

* Fixed adding summary columns not taking columns into account.
* Added old summary column format support.
* Fixed support for bigint type.
* Fixed exception thrown onHit.

## 1.2.0 (2023-11-02)

### Features

* [#2208](https://github.com/datagrok-ai/public/issues/2208): Implemented SmartForms.
* Added image rendering on double clcik in ImageUrl.

### Bug Fixes

* Removed URL rendering in ImageUrl.
* Fixed sparklines visualizations and made them configurable.

## 1.1.33 (2023-08-04)

### Bug Fixes

* Fixed scroll on pinned columns.
* [#2117](https://github.com/datagrok-ai/public/issues/2117): Pinned columns are not displayed properly with some row source options for the table (MouseOverRow, FilteredSelected, Selected, SelectedOrCurrent).

## 1.1.32 (2023-08-02)

### Bug Fixes

* Added tests for pinned columns

## 1.1.31 (2023-08-01)

### Bug Fixes

* [#2117](https://github.com/datagrok-ai/public/issues/2117): Pinned columns are not displayed properly with some row source options for the table

## 1.1.30 (2023-08-01)

## 1.1.29 (2023-07-31)

### Bug Fixes

* Fixed mouse wheel handler for pinned columns

## 1.1.28 (2023-07-28)

### Bug Fixes

* Fixed undefined is not iterable problem

## 1.1.27 (2023-07-17)

### Features

* Added structure editing support for pinned columns.

### Bug Fixes

* [#1912](https://github.com/datagrok-ai/public/issues/1912): Can't edit values in pinned column.
* GROK-13348: Can't edit data via a single click for pinned boolean column.
