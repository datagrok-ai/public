# Statistics — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `stats_viewer/CLAUDE.md` "Automation
surface" and `statistics-ui.md` / `statistics.md`. Back Color is in `statistics-back-color.feature`.

## Still open

- Reordering stat columns (dragging the `max` header past `min`): needs the inner grid's `column order` re-exported in
  statistic names.
- Resizing stat columns: needs `column resizer <stat>` areas and `column width of <stat>` readings re-exported.
- Full screen with Alt+F and back: needs a `full screen` reading on the base viewer status.

## Won't do

- A Histograms column dropped when its column is removed, a renamed column, a calculated column added: the rules every
  viewer shares (`../MISSING.md`), xamgle `d4 | Viewers | Statistics rebuilds once for a burst of removed columns`.
- Statistics over no rows, a constant or an all-empty column (core#1125): ddt `stats_test.dart` and `med_test.dart`.
