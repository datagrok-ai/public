---
feature: heatmap
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.heat-map, viewers.filters.histogram]
realized_as:
  - heatmap-spec.ts
related_bugs:
  - id: GROK-20619
    status: open
  - id: GROK-11515
    status: fixed
---

# Heat map (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Heat map** from **Toolbox > Viewers**

## Add the viewer

1. Click the **Heat map** icon in **Toolbox > Viewers**
2. The map is drawn and **Is Heatmap** is on

## Heatmap colours

1. Uncheck **Heatmap Colors** in **Context Panel > Misc** — the cells should stop
   being colour-filled. They do not change today (GROK-20619, guarded by
   `knownOpenBug`). Check it again.
2. Flip **Global Color Scaling** — the cells are recoloured, then restored

## Column labels

1. Set **Col Labels Orientation** (Context Panel > Style) to *Vert*, then *Horz*,
   then back to *Auto*

## Columns and scrollbars

1. Set **Max Heatmap Columns** to *3* — fewer columns are drawn; restore *100*
2. Uncheck **Show Heatmap Scrollbars** — the range sliders disappear; check it again

## Grid mode

1. Uncheck **Is Heatmap** — the viewer redraws as a plain grid
2. Check it again — the heat map comes back

## Row height

1. **Row Height** should not be offered here at all: it is a grid-only option
   (its tooltip says so) and the drawing does not follow it in heatmap mode. The
   property is still shown today (GROK-20619, guarded by `knownOpenBug`); once it
   is hidden, this check goes loud and the step can be dropped.

## Interaction

1. Click a cell — its row becomes current
2. **Alt+drag** over an area — the view zooms in and the vertical slider's window
   shrinks
3. Double-click the vertical slider — the view resets to the full extent

## Filtering

1. Filter the table to `SEX = M` — the heat map redraws; reset the filter

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`heatmap-spec.ts`. Kept verbatim so no scenario is lost.

### Table switching

> Manual

> Note: requires spgi-100 dataset (open twice to get two tables for switching).
> Setup: Close all, open spgi-100.csv twice, go to the first table view, add Heat map, open Context Panel.

1. Set the Table property to the second spgi-100 table — heat map re-renders against it
2. Set the Table property back to the first spgi-100 table

### Range slider navigation

> Manual

1. Drag the horizontal range slider to the left — content scrolls left
2. Drag the horizontal range slider to the right — content scrolls right
3. Double-click the horizontal range slider — resets to full range
4. Drag the vertical range slider up — content scrolls up
5. Drag the vertical range slider down — content scrolls down
6. Double-click the vertical range slider — resets to full range

### Column sorting

> Manual

1. Double-click the AGE column header — rows sort ascending by AGE
2. Double-click the AGE column header again — rows sort descending
3. Double-click the AGE column header once more — sorting resets

### Selection interaction

> Manual

1. Click a row in the heat map — it becomes the current row
2. Shift+drag to select a block of cells — selected rows are highlighted
3. Press Esc — selection clears

### Color scheme customization

> Manual

1. Open properties, locate Linear Color Scheme
2. Change the linear color scheme to a different option
3. Open properties, locate Categorical Color Scheme
4. Change the categorical color scheme to a different option

### Draw every row

> Manual

1. Open properties, verify Draw Every Row is unchecked
2. Check Draw Every Row
3. Scroll through the heat map
4. Uncheck Draw Every Row

### Layout save and restore

> Manual

1. Set Col Labels Orientation to Vert
2. Check Global Color Scaling
3. Save the layout
4. Reset Col Labels Orientation to Auto and uncheck Global Color Scaling
5. Apply the saved layout — orientation and scaling restore
6. Delete the saved layout

### Layout saving with Is Heatmap toggle

> Manual

> Note: uses spgi-100 dataset.
> Setup: open spgi-100.csv, add Heat map.

1. In the Property Pane, set Max Heatmap Columns to 200
2. Uncheck Is Heatmap — viewer switches to grid mode
3. Save the layout
4. Reset Max Heatmap Columns to 100 and check Is Heatmap
5. Apply the saved layout
6. Check Is Heatmap — re-enable heatmap mode
7. Delete the saved layout

### Right-click content panning

> Manual — see `heatmap-ui.md`

1. Zoom in using Alt+drag, then right-click and drag on the content area — content pans
2. Release right-click — panning stops, content stays at new position

---
{
  "order": 14,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
