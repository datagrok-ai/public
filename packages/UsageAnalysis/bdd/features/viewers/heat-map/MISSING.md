# Heat map — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against the grid's "Automation
surface — heat map mode" (`core/client/d4/lib/src/viewers/grid/CLAUDE.md`), the old TestTrack checklist
(`files/TestTrack/Viewers/Heatmap/heatmap.md`) and the library's rules; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- Which tab a stacked heat map comes back on (lead, 2026-10-19): a layout improvement, core#931;
  `heat-map-layout.feature` switches to the Heat map tab itself until it lands.
- Max Heatmap Columns counts only the columns the user shows, and a hidden column stays hidden when the cap rises
  (core#1202): Dart `Viewers | Grid | Max heatmap columns counts only the columns the user shows`.

## Still open

- **Linear / Categorical Color Scheme** (Style): no reading says which scheme a band is drawn with; needs a
  `linear color scheme` / `categorical color scheme` reading (the scheme's colours, `#rrggbb` joined).
- **Draw Every Row**: nothing reports which path the last frame took; needs a `draws every row` reading (the
  `renderEveryRow` the frame used, `grid_core.dart`).
- **The tooltip over a heat map band (public#732)**: the dense mode reports one band per column, not the row under the
  pointer; needs a `mouse-over row` reading in heat map mode, or a zoom deep enough for per-cell areas, walked by hand
  first.
- **`column resizer <col>` in heat map mode**: `grid_status.dart` reports the band while `ColumnResizing.check` refuses
  heat map mode (GROK-18672); it should be reported only where `check` accepts it.
- **Right-click panning** (`heatmap-ui.md`): the `row order` reading moves with a pan, but no step drags with the
  right mouse button; needs a library gesture.
- **The viewer's error in the grid status**: with Max Heatmap Columns 0 the map shows "No columns to show: Max Heatmap
  Columns is 0" (core#928), yet `getWidgetStatus().error` is null; `grid_status.dart` should report it, then the hint is
  claimed with `heat map viewer should report the error "…"`.

## Won't do

- **What a layout stores** (Col Labels Orientation, Global Color Scaling, the scroll window, public#2926): an API test;
  `heat-map-layout.feature`'s public#2927 scenario reads the restored window as a person sees it.
- **Sorting the spreadsheet next to the heat map**: each grid holds its own order; not a heat map claim.
- An Alt-click zooming to one row: incidental. An Alt-drag on a sorted map, Sort > Custom with the empty value first:
  unmeasured, and the custom order dialog is the grid's. AGE removed under the map: the shared viewer rule (Dart
  `Column removal`).
