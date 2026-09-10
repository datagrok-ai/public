# GIS changelog

## v.next

* Map: Fixed `isRenderPending` never clearing under Render Type = heatmap - the hidden WebGL marker layer keeps `renderer.ready === false` and OpenLayers never renders it again, so every settle waited out its cap
* Map: Added the automation surface — `getWidgetStatus` reports `view`, `layers panel`, `layer "<name>"`, `layer "<name>" visibility`, `zoom in` / `zoom out` and `point <row>` hit areas plus `layers`, `visible layers`, `layer "<name>" visible`, `zoom`, `centre`, `markers`, `rows shown` and `render type`; added `isRenderPending` / `onRendered`, settled on the marker layer rather than on base-map tiles
* Fixed crash on closing Map viewer caused by pointer events firing after disposal
* Layers panel: Fixed page freeze when toggling layer visibility (e.g. Markers GL) — caused by exponential re-subscription of `onCurrentCellChanged` on every layers list refresh
* KML/KMZ: Fixed `SyntaxError: Unexpected end of JSON input` in `GisAreaGridCellRenderer` when the gisObject cell is empty (style-only KML rows or unsupported geometries)

## 1.3.1 (2025-03-30)

* Fixed build issues

## 1.2.1 (2024-03-20)

* Remove JSON handling leaving only geo/topo json. 

## 1.1.5 (2024-03-01)

### Bug Fixes

* Fixed 'Cannot read properties of null (reading 'hasOwnProperty')' error when closing a viewer.

## 1.1.4 (2023-07-24)

### Bug Fixes

* Fixed the test dataset location in the package.

## 1.1.3 (2023-06-13)

### Breaking change

* Removed zipcode detector since it's not used.
