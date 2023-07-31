# Charts changelog

## 1.1.0 (WIP)

This release focuses on adding new functionality and improving the existing one.

*Dependency: datagrok-api >= 1.16.0*

### Features

* Switched to class decorators for registering viewers
* [#1974](https://github.com/datagrok-ai/public/issues/1974): Timelines: Set the default date format based on the time difference between event start and end values
* [#2096](https://github.com/datagrok-ai/public/issues/2096): Sunburst: Implemented selection capability on Ctrl + Click

### Bug Fixes

* [#1968](https://github.com/datagrok-ai/public/issues/1968): RadarViewer does not allow changing columns in SpiderChart:
  * Show actual number of columns that are selected
  * RadarViewer contains non numerical columns
* [#1962](https://github.com/datagrok-ai/public/issues/1962): RadarViewer errors on bigint datatype
* Surface plot fixes (format, axis labels)

## 1.0.24 (2023-05-26)

This release focuses on stability.

*Dependency: datagarok-api >= 1.15.0*

### Bug Fixes

* [#1973](https://github.com/datagrok-ai/public/issues/1973): Timelines: Visualization is not shown if the Split by column contains long category names
