# GIS

Geospatial functionality for Datagrok: geocoding, GIS semantic types, and the **Map** viewer —
`GisViewer` (`src/gis-viewer.ts`), a `DG.JsViewer` wrapping OpenLayers 6.15.1 (`src/gis-openlayer.ts`)
that draws a table's latitude/longitude columns as WebGL markers or a heat map over a base map.
Registered as `Map` (`src/package.ts:441`), so a feature reaches it as **`map viewer`**.
Note that the Dart `Shape Map` viewer is a different thing entirely — an SVG choropleth, no
OpenLayers, no tiles.

| File | Purpose |
|---|---|
| `src/gis-viewer.ts` | The viewer: properties, coordinate extraction, render modes, the three automation signals |
| `src/gis-openlayer.ts` | The OpenLayers wrapper: the map, its layers, sources, styles and interactions |
| `src/gis-mapcontrols.ts` | The two map controls: the layers button and the layers panel |
| `src/gis-viewer-status.ts` | The areas and readings `getWidgetStatus` reports |

## A pixel claim on this viewer is not a product claim

`initMap` builds two tile layers before anything of ours is drawn (`gis-openlayer.ts:356-357`): a
`BingMaps` aerial layer (`:806-824`, with an API key in the source) and an `OSM()` base layer
(`:825-843`). **Both fetch from the public internet.** Any assertion that hashes a screenshot of
this viewer, or counts its ink, is an assertion about a tile server's availability and latency, not
about Datagrok. Offline it fails; online it is non-deterministic, because a tile landing a frame
later changes the picture.

So features against this viewer read state — `layer "…" visible`, `zoom`, `centre`, `markers`,
`rows shown` — and point at named areas. They never compare pictures. There is a further reason the
picture is not even well defined here: the tile layers and the WebGL marker layer render onto
**separate canvases** inside `.ol-viewport`, so "the canvas" of this viewer is whichever one comes
first in the DOM.

## The layers panel is the OL control, not `divLayersList`

`GisViewer` builds a `panelLeft` box holding a `divLayersList` (`gis-viewer.ts:60, 344-345`), but
`initUi` appends only `viewerContainer` to the root (`:360`) — the split layout that would have held
`panelLeft` is commented out (`:355-358`). `panelLeft`, `panelTop`, `panelBottom` and `divLayersList`
are therefore **never in the document**, and `updateLayersList` returns immediately anyway while
`isShortUI` is true (`:420-421`), which is its only mode.

The layers panel the user sees is `PanelLayersControl` (`gis-mapcontrols.ts:12-206`), an OpenLayers
control added in `initMap` (`gis-openlayer.ts:344-346`) whose element is `.panel-layers`. Inside it
is a real `DG.Grid` over a one-row-per-layer dataframe with columns `vis`, `name`, `exp`, `del`
(`gis-openlayer.ts:681-708`), row height 24. The panel starts hidden (`gis-openlayer.ts:385`) and the
layers button toggles it.

That grid is where the layer hit areas come from: the status asks it for its own cell bounds, so the
row a layer sits in is the map's layer order, and there is no 24-pixel arithmetic anywhere.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, `canvas` (the first canvas inside `.ol-viewport`), and `layers panel` — the latter
only while the panel is showing. Hit areas are in CSS px of that canvas.

| Area | What it is, and when it is absent |
|---|---|
| `view` | The `.ol-viewport` box — the map itself, excluding nothing, since the controls float over it. |
| `zoom in` / `zoom out` | The default `ol-zoom` control's two buttons. |
| `layers panel` | `.panel-layers`. **Absent while the panel is hidden** (its `visibility`), which is how it starts. |
| `layer "<name>"` | The layer's row in the panel grid, from the `vis` cell's left edge to the `name` cell's right edge. Absent while the panel is hidden, and for a layer whose `layerName` is unset. Two layers with the same name collapse to one area — the last one wins. |
| `layer "<name>" visibility` | Just the `vis` checkbox cell of that row — the control that toggles the layer. This is the one to click; the middle of `layer "<name>"` lands on the name, which only makes the layer current. |
| `point <row>` | A box of `2 × markerDefaultSize` px centred on where the map projects that row's coordinate, for at most the first 20 filtered features that fall inside the viewport. `<row>` is 1-based. Absent for a row outside the viewport, for a row with no point geometry, and for all rows before the map's first frame (`getPixelFromCoordinate` has no frame state to work from). The box is a pointing target of the default marker size, not the marker's drawn size — a size column changes what is painted, not this. |

| Reading | Reads from |
|---|---|
| `layers` | `ol.getLayersList().length` — every layer on the map, base maps included. |
| `visible layers` | How many of them are visible. |
| `layer "<name>" visible` | That layer's `getVisible()`. Replaces walking the map's layer collection from the test. |
| `zoom` | `getView().getZoom()`, rounded to 2 decimals. Fractional during a zoom animation, so compare after a settle. |
| `centre` | `lon, lat` of `getView().getCenter()` in EPSG:4326, 3 decimals — the map's projection is EPSG:3857 (`gis-openlayer.ts:347-352`), so the raw centre is in metres and is not what a feature should read. |
| `markers` | `features.length` — the features currently handed to the map, which is the filtered set. |
| `rows shown` | `dataFrame.filter.trueCount`. Equal to `markers` in the normal case; they **diverge** when the latitude or longitude column is unset, because `getCoordinates` clears the feature arrays and returns (`gis-viewer.ts:775-779`) while the filter is untouched. That divergence is the honest witness for "geo columns detected". |
| `render type` | `markers`, `heatmap` or `both`. |

`shortcuts` and `events` are empty, and `error` is always `null`: the viewer reports an init failure
through `grok.shell.error` (`gis-viewer.ts:406`), a balloon, and has no error surface of its own.

## `isRenderPending` and `onRendered`, and why not `rendercomplete`

The obvious signal is the map's `rendercomplete`, and it is the wrong one. OpenLayers computes it as
`!tileQueue.getTilesLoading() && !tileQueue.getCount() && !getLoadingOrNotReady()`
(`ol/PluggableMap.js:1319-1326`), so it waits on **every tile in the viewport**. On a stand with
internet that turns every settle into a tile download; on a stand without one it depends on tiles
failing promptly. (They do drain the queue — `TileQueue.handleTileChange` counts `ERROR` as done,
`ol/TileQueue.js:88-98` — and the Bing metadata request, which never resolves offline because
`BingMaps` passes no errback to `jsonp` (`ol/source/BingMaps.js:178`, `ol/net.js:31-57`), does not
block it either: `getLoadingOrNotReady` reads `source.loading`, a flag tile sources never set
(`ol/source/Source.js:92`). So the event is not certain to hang — it is certain to be slow and to
depend on the network, which is enough to disqualify it.)

Instead the viewer counts a render as done on the map's next `postrender` — dispatched every frame
regardless of tiles (`ol/PluggableMap.js:1318`) — on which the WebGL marker layer's renderer reports
ready (`gis-viewer.ts`, `renderFinished`). The marker layer is the product; the base map is not.

The flag goes up when a change arrives, not when the debounced handler runs: the raw
`selection.onChanged` and `filter.onChanged` streams raise it ahead of their 100 / 200 ms debounces,
and `render` raises it too. `render` then asks the map for a frame in its `finally`, so the
`postrender` that clears the flag is always a frame requested after the work — a map that would
otherwise stay idle still settles.

## Writing features against this viewer

**Open the layers panel before asking for a layer area.** It starts hidden; `layers panel` and every
`layer "…"` area are absent until the layers button is clicked. The `layer "…" visible` **readings**
need no panel at all — prefer them.

**`layers` counts the base maps.** On a fresh Map viewer over a geo table the list is Bing sat,
BaseLayer, Markers GL, Markers GL Selection and Heatmap — asserting a bare number will break the
first time a layer is added. Assert the layer you care about by name.

**Fixture: `earthquakes`** (`System:DemoFiles/geo/earthquakes.csv`, registered as the bdd dataset
`earthquakes` in `libraries/bdd/bindings/platform/datasets.ts:10`) — 2426 rows; DateTime, Latitude, Longitude, Depth, Magnitude, MagType and more.
MagType is Mw 2312, Mb 72, Me 25, ML 7, Unk 6, Ms 4.

**Do not write a `@needs-network` picture claim.** If a scenario cannot be made without comparing
tiles, it does not belong in the feature.

## Rebuilding after a change here

GIS is a published package: `npm run build` (`grok api`, `grok check --soft`, webpack) then
`grok publish <alias>`, or `npm run build-gis` for the webpack step alone. `rxjs` is a webpack
external, so the platform provides it at runtime. Nothing in the client or the Dart core has to be
rebuilt for these keys.
