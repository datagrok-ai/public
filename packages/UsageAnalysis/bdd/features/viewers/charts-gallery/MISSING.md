# Charts gallery (Surface plot, Multiplot) — what the features cannot say yet

Written during the Charts gap round of 23 Oct 2026, against `public/packages/Charts/src/viewers/`
and the library's rules (`public/libraries/bdd/CLAUDE.md`); trimmed in the BDD PR round of 8 Nov 2026.

## Still open

- **The Surface plot has no status** (`surface-plot.ts`, echarts-gl): its Projection and Wireframe claims are the
  properties read back. Needs `getWidgetStatus` readings (`x column`, `y column`, `z column`, `points` after
  `plotFilter`, `projection`, `wireframe` of the applied option) and `isRenderPending` until echarts' `finished`.
- **The Multiplot has no feature and no status**: needs a status listing its inner viewers and exposing them as inner
  widgets.

## Won't do

- The Globe under a filter-panel category and Color By picked in the context panel (PR 4516): the points following
  the table's filter are claimed in `charts-gallery.feature` through the filter itself (the panel's click is generic);
  Color By read back is an echo.
