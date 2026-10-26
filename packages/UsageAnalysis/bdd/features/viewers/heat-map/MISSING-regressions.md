# Heat map — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

- **The grid status does not publish the viewer's error** (measured 26 Oct 2026 on localhost, master).
  With Max Heatmap Columns 0 the map shows "No columns to show: Max Heatmap Columns is 0" in its
  `.d4-viewer-error` (core#928, `grid_core.dart` sets it as the viewer error), yet
  `getWidgetStatus().error` is null, so `heat map viewer should report the error "…"` cannot see it and
  `should report no error` passes over it. `heat-map-regressions.feature` claims the hint as text shown
  in the viewer (`heat map viewer should show the text …`). Needs `grid_status.dart` to report the
  viewer's error, as the other viewers' statuses do; then the hint can be claimed as:

  ```gherkin
  When user sets "maxHeatmapColumns" property of heat map viewer to "0"
  Then heat map viewer should report the error "No columns to show: Max Heatmap Columns is 0"
  When user sets "maxHeatmapColumns" property of heat map viewer to "100"
  Then heat map viewer should report no error
  ```

  `heat-map-columns-and-sorting.feature`'s GROK-17776 scenario asserts `should report no error` right
  after Max Heatmap Columns 0; it passes only because of this gap, and fails once the status reports
  the error.
