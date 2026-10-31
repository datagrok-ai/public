# Heat map — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against
the grid's "Automation surface — heat map mode" (`core/client/d4/lib/src/viewers/grid/CLAUDE.md`), the
old TestTrack checklist (`files/TestTrack/Viewers/Heatmap/heatmap.md`, its "Manual" sections) and the
library's rules (`public/libraries/bdd/CLAUDE.md`). The three new features hold what could be said:
`heat-map-layout.feature` (public#2927, public#2901, GROK-16239), `heat-map-columns-and-sorting.feature`
(the Columns property and its dialog, core#376; header sorting, public#182; no column resize in heat map
mode, GROK-18672; Max Heatmap Columns 0, GROK-17776) and `heat-map-interaction.feature` (the scroll
slider handles, a Shift-drag selection and Escape, the Table property).

## Still open

- **Linear / Categorical Color Scheme.** `GridLook.linearColorScheme` and `categoricalColorScheme`
  (Style) have no reading: the status reports `heatmap colors` and `global color scaling`, not which
  scheme a band is drawn with. A repaint of the `column AGE` band after a scheme change would pass for
  any change at all, and the API property step cannot write the list (`"[4278255360, 4278255360]"`
  became the two strings `"[4278255360"`, `"4278255360]"`); the UI path is the colour-scheme picker in
  the property grid, which no step drives. Needs a `linear color scheme` / `categorical color scheme`
  reading (the scheme's colours, `#rrggbb` joined), and then a claim that a band contains the scheme's
  end colour (`the "column AGE" area … should contain the color …`). The categorical scheme applies
  only to columns of 100+ categories (`grid_look.dart`), which demog-1000 has none of: USUBJID is a
  1000-category string column, so the claim would go on its band.
- **Draw Every Row.** `drawEveryRow` decides whether a sub-pixel heat map blends every row or samples
  one per pixel (`grid_core.dart`, `renderEveryRow`); nothing reports which path the last frame took,
  and on demog-1000 the two pictures differ only in blending, so "repainted" is the only claim — one
  any repaint satisfies. Needs a `draws every row` reading (the `renderEveryRow` the frame used) so the
  property can be paired with what the frame did; the md's "scroll through the heat map" adds nothing a
  reading would not.
- **The tooltip over a heat map band (public#732).** public#732 asked for a tooltip that shows the full
  USUBJID on a zoomed heat map. A hover at the centre of the `column USUBJID` band
  (`user hovers over the "column USUBJID" area of heat map viewer`) showed no visible `.d4-tooltip`
  on the stand, and the dense mode reports one band per column, not the row under the pointer, so the
  tooltip steps have nothing to compare with. Needs either a `mouse-over row` reading in heat map mode
  or a zoom deep enough for per-cell areas (row height over 8 px) before the hover; not written until
  the person's path (zoom in, hover a USUBJID cell) is walked by hand.
- **`column resizer <col>` in heat map mode.** `grid_status.dart:62` reports the resizer band for every
  header whenever `showColumnLabels` is on, but `ColumnResizing.check` refuses heat map mode
  (`!grid.look.isHeatmap`, GROK-18672), so the status names a drag target that does nothing.
  `heat-map-columns-and-sorting.feature` drags it on purpose (the width stays 30) and pairs it with the
  same drag in grid mode; still, the area should be reported only where `check` would accept it, as
  the grid's own CLAUDE.md promises ("the same band `_checkResizeHandleAt` accepts").

- **Right-click panning** (`heatmap-ui.md`): the reading exists — `row order` starts at the window's
  first row (`grid.minRow`), so a vertical pan of a zoomed map moves it while `y scroll span` stays —
  but no step drags with the right mouse button (every drag step presses the left one; `{key}` holds
  keyboard modifiers only). Needs a library gesture, `user drags the {string} area of {widget} by
  {int} pixels to the {word} with the right button`; not written as a package step because it is not
  a composition of existing calls.

## Stale lines in `core/client/d4/lib/src/viewers/grid/CLAUDE.md`

- **Lines ~363-366** ("Max Heatmap Columns does **not** move this reading on a live viewer … leaves all
  11 columns shown"): stale since the look binds its grid on every viewer (grid_core.dart, 2026-09-21,
  core#157). `heat-map.feature` and `heat-map-columns-and-sorting.feature` set the property on a live
  viewer and read `columns shown` move to 3 (or 0) and back to 11.
- **Lines ~389-391** ("the property row exists and is `aria-disabled` — read that state, do not assert
  the row's absence"): stale since GROK-20619 — a property gated on the non-editable `isGrid` is left
  out of the panel (`prop_grid.dart` `isLockedOff`), and `heat-map.feature` asserts the Row Height row
  `absent`.
- **Line ~326** (`column resizer <col>` described as the band `_checkResizeHandleAt` accepts): true in
  grid mode only, see above.

## Answered

- **Which tab a stacked heat map comes back on** (lead, 2026-10-19: not a heat map bug — an improvement for
  layouts, core#931). After the layout loads, the panel shows the Histogram tab, not the Heat map tab that
  was active at the save; a layout does not keep the active tab of stacked viewers. `heat-map-layout.feature`
  switches to the Heat map tab itself before claiming it; once core#931 lands that switch can go.

## Won't do

- **What a layout stores** (Col Labels Orientation and Global Color Scaling restored from a layout, the
  md's "Layout save and restore"; the scroll window restored, public#2926): by the lead's decision a
  layout scenario is a regression of visible breakage; the round trip of a setting is an API test.
  `heat-map-layout.feature`'s public#2927 scenario reads the restored window only as what a person
  sees (as many columns on screen as before, not one).
- **Sorting the spreadsheet next to the heat map**: a double-click on the table view's grid header
  sorts that grid only — the heat map keeps the table's order (`row order` stays 1..10), since each
  grid holds its own order. That is how every grid behaves, not a heat map claim.

## Closed issues the regression round of 26 Oct 2026 could not claim

- **The grid status does not publish the viewer's error** (measured 26 Oct 2026 on localhost, master).
  With Max Heatmap Columns 0 the map shows "No columns to show: Max Heatmap Columns is 0" in its
  `.d4-viewer-error` (core#928, `grid_core.dart` sets it as the viewer error), yet
  `getWidgetStatus().error` is null, so `heat map viewer should report the error "…"` cannot see it and
  `should report no error` passes over it. `heat-map-columns-and-sorting.feature` (GROK-17776, core#928) claims the hint as text shown
  in the viewer (`heat map viewer should show the text …`). Needs `grid_status.dart` to report the
  viewer's error, as the other viewers' statuses do; then the hint can be claimed as:

  ```gherkin
  When user sets "maxHeatmapColumns" property of heat map viewer to "0"
  Then heat map viewer should report the error "No columns to show: Max Heatmap Columns is 0"
  When user sets "maxHeatmapColumns" property of heat map viewer to "100"
  Then heat map viewer should report no error
  ```
