---
feature: map
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - map-spec.ts
related_bugs: []
---

# Map viewer (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/geo/earthquakes.csv**
3. Add the **Map** viewer from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, the Context Panel
property grid, the map's own layers panel and zoom buttons, and real mouse and
keyboard input on the map. The JS API is used only to read back state that the UI
does not expose as text (OpenLayers layer visibility, zoom level, selected-row count).

## Add the viewer and auto-detect the geo columns

1. Click the **Map** icon in **Toolbox > Viewers** — a map with points appears
2. On the **Context Panel > Data**, check **Latitude** shows *Latitude* and
   **Longitude** shows *Longitude* — the columns are detected automatically
3. Check the layer stack: **Markers GL** is visible, **Heatmap** is hidden

## Color and size coding

1. On the **Context Panel > Data**, set **Color** to *Magnitude* — the panel shows
   *Magnitude* and the map repaints
2. Set **Size** to *Depth* — the panel shows *Depth* and the map repaints

## Marker sizing

1. Go to **Context Panel > Markers**
2. Set **Marker Min Size** to *12* — the markers grow, the map repaints
3. Set it back to *2*

## Layers panel

1. On the map, click the **Map layers** icon — the layers panel opens
2. Check **Heatmap** — the heatmap layer becomes visible and the map repaints
3. Uncheck **Markers GL** — the markers disappear, and the UI stays responsive
   (this used to freeze the page)
4. Restore both layers to their original state and close the panel

## Render type

1. Go to **Context Panel > Misc**
2. Set **Render Type** to *heatmap*, then *both*, then back to *markers*
3. Each mode paints a visibly different map

## Zoom controls

1. Click **+** on the map — the zoom level grows by exactly 1 and the map repaints
2. Click **−** — the zoom level returns to where it started

## Selection

1. **Ctrl + drag** a rectangle over the map — the enclosed points are selected,
   the selection is reflected in the table, and the map repaints
2. Press **Escape** — the selection clears

## Filtering

1. Filter the table down to `MagType = Mw`
2. Fewer rows than the full table remain and the map repaints — the map shows
   only filtered rows
3. Reset the filter

## Tooltip

1. Go to **Context Panel > Misc** and enable **Show Tooltip**
2. Hover a marker — a tooltip appears listing the row's column values
   (including *Latitude* and *Magnitude*)

## Closing the viewer

1. Move the pointer over the map, then click **Close** on the viewer title bar
2. The viewer is gone and no error is raised — pointer events arriving after
   disposal used to crash the viewer

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`map-spec.ts`. Kept verbatim so no scenario is lost.

### Heatmap layer settings

> Manual

11. Go to **Property Pane > Heatmap** and use sliders to change the **Heatmap Radius** and **Heatmap Blur** properties**.**

### Mouse navigation

> Manual

14. Go to the map viewer and use the mouse wheel to zoom in and out.
16. Use **Mouse Drag** on the map to move it.

### Selection from the grid

> Manual

19. Go to the grid and select some rows. Check the selection on the map.
20. Go to **Property Pane > Markers** and change the **Selected Color** property.

### Custom geo layers

> Manual

22. Drag and drop ***doc.kml*** right to the map. A **doc.kml** layer should appear.
23. Go to the viewer and uncheck the **doc.kml** checkbox.
24. Delete the **doc.kml** layer.
26. Drag and drop the ***ecoregions2.geojson*** file (repeat for ***boothbay2018.kmz, doc.kml, uk.topojson***). The file with a map viewer should appear.
27. Click objects on the viewer, and check the map viewer interaction.

---
{
  "order": 11,
  "datasets": ["System:DemoFiles/geo/earthquakes.csv"]
}
