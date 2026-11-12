# Curves features — what is missing

## Resolved

- **Multi curve viewer readings** (2026-10-18) and **its axes** (GROK-20661): `getWidgetStatus()` reports the curves,
  the legend and `x range` / `y range` / `points outside the axes`; `multi-curve-viewer` and `curves-legend` claim them.
- **core#1246** (Data to Curves missing from the table view's Data menu): fixed in core (`table_view.dart` `fillMenu`
  takes functions with `requiresTableView`); `data-to-curves.feature` picks it from the menu.

## Still open

Nothing planned.

## Won't do

- **Curve cell readings in the grid and in Context Panel > Chart, points as hit areas in a cell**: what a renderer
  draws, the package tests (legend-tests, panel-renderer-tests, converter-tests); ruling of 2026-10-06.
- **public#4294** (a new Multi curve viewer's Show Legend and legend): the `curves-legend` guide claims both.
- **Data to Curves**: what each Fitted Curves row holds (group, run, maximum inhibition) is the package test 'Data to
  Curves'; the other prefilled roles are one name rule; the Parent Table side is the converter's.
