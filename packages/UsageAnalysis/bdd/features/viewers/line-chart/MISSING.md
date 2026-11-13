# Line chart — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `line_chart/CLAUDE.md` "Automation
surface"; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- A chart type no renderer draws (`Bar Chart`, an unknown name from a layout) is drawn as a line instead of throwing
  on every frame: Dart `Line chart draws a chart type it has no renderer for as a line`.

## Still open

- **Show Mouse-over Category** (`line-chart.md`): the property works (walked by hand); needs a `mouse-over category`
  reading in `line_chart_status.dart`.
- **A click on a line selects its category; Control adds, Control+Shift removes (public#2608, public#2360)**: needs a
  `line <category>` area (a point on a drawn segment that `hitTestLineCategory` accepts).
- **A split or a Y column removed on the plot shows in the property panel at once (public#2951, public#2875,
  public#2623, public#2880)**: needs `user removes the {string} split of {widget}` and `… Y column of {widget}` steps (the
  cross of that column's chip).
- **A logarithmic axis over non-positive values shows its warning (public#3701)**: needs `x warning` / `y warning`
  areas, as the density plot reports.
- **After a category filter the X labels stay over their points (public#1456)**: needs `x label <category>` areas and a
  geometry step that centres one area over another.
- **A custom category order is what the axis and the points follow (public#1551)**: needs a step that sets a column's
  custom category order, and the `x label <category>` areas.
- **Multi Axis with Y Global Scale puts a constant series on its own tick (core#823)**: needs a `y of point of row <n> of
  "<col>"` reading.
- **Resizing does not join the first and last points (public#1671); the SPC average label keeps off the axis labels
  (public#4013)**: need `segments of "<col>"` and an `average label` area.

## Won't do

- Every row filtered out under a drag, the Y column removed (core#1100): the shared viewer rules — Dart `Data
  Manipulation | Line chart`, `Line chart with a Y column missing from its table hit-tests and maps points`.
- The regression line after a table swap, a stacked area under a filter, the SPC average of the filtered rows, the
  aggregated markers, the restore timer of a closed chart (core#1100): the Dart line chart regTests.
- X Min / X Max outside the data, Line Width 0 or −3, a logarithmic Y over non-positive values, whiskers over one-row
  groups, ±Infinity, all-null columns, Min above Max, a layout whose split column is gone: not UI checks, or the shared
  viewer rules.
