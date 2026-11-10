# PC plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `pc_plot/CLAUDE.md` "Automation surface"
and `pc-plot-ui.md`. What could be said is in `pc-plot-hover-and-layout.feature`.

## Still open

- Show Mouse Over Row Group ("Mouse-over and row group", item 2): `lines drawn` counts only the plot's own lines; needs
  a `mouse-over group lines drawn` reading (the overlay's group pass in `pc_plot_core.dart`).
- Undock and dock again ("Viewer basics", items 1-2): no step floats a viewer out of the dock manager.
- General > Save to Gallery: it leaves a server entity, which needs a gallery entry swept at feature end.
- core#1339 — General > Save as PNG logs three 404s for `/font/MaterialIcons-Regular.*` (`google_fonts.css`); claim when
  fixed: `And no errors should have been logged` after the download in `pc-plot.feature`.

## Won't do

- A grid column header dropped on the plot: the plot takes no foreign drop (`pc_plot_core.dart:125`, its own axis labels
  only); adding an axis that way is a feature request.
- An Axis Font the canvas cannot parse (core#967): reachable only through the API or a layout; the panel's font editor
  makes valid fonts.
- A layout with Line Width 0 and density on, hovered (core#967): ApiTests `AI: Viewers: PcPlot Extras`.
- A Log Columns axis over zeros or negatives (core#1203): xamgle `d4 | Viewers | Fixed issues | PC plot Log Columns with
  zero or negative values`.
- Bins, Horz Margin and the line widths past their bounds, an axis column removed or renamed, a layout naming a missing
  column, no rows and deleted rows: the rules every viewer shares (`../MISSING.md`).
- One axis left (the "at least 2 columns" message), a constant column under Normalized with a density: low value.
