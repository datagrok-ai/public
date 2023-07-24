# Curves changelog

## 1.1.0 (2023-07-21)

This release focuses on improving analysis stability and usability.

*Dependency: datagarok-api >= 1.16.0*

### Features

* Improved curve fitting and confidence intervals:
  * Added a semantic type detector for "fit" to render fitted curves.
  * Fixed performance issues during fitting.
  * Automated data sorting during fitting.
* Implemented interactive parameter recalculating.
* Added different marker types.
* Improved axes cell rendering.
* Added the possibility to run curves by user-defined Javascript function.
* [#2100](https://github.com/datagrok-ai/public/issues/2100): Made cell renderer edit mode a resizeable window.
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improved curves properties and rendering: 
  * Don't render curves before the axes start.
  * Changed margins in grid cell for axes
  * Added the capability to control parameter bounds.
  * Fixed curves rendering by inverting the y while converting the coordinates.
* [#2105](https://github.com/datagrok-ai/public/issues/2105): Outliers are now shown in red color and with a much bigger size.
