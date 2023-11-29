# PowerGrid changelog

## 1.2.5 (WIP)

### Features

* Added default bottom position for Forms viewer.

## 1.2.4 (2023-11-27)

### Features

* Added Forms viewer icon.

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
