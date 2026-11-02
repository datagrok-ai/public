# Pie chart edge cases — what the features cannot say yet

The edge cases of this viewer that the features could not state with the current step library and the expected behaviour documented today; the claimed ones are `@edge-cases` scenarios of their own, after the master scenarios of `pie-chart-filter-and-auto-legend.feature`, and the known failures in `pie-chart-known-failures.feature`. No library, binding or code change was made for them.

## Still open

**Decided** (the rules in `../MISSING.md`)

1. A renamed category column is followed, and every row filtered out leaves no slice: requested below as API tests.
2. Max Radius typed in the context panel above its maximum of 500 (600) is not taken; Start Angle typed as −30 or 400 is the same. Both are known failures of core#1200, the scenario outline of `pie-chart-known-failures.feature` (one Examples row per value; back into `pie-chart.feature` once fixed).
3. The category column removed: no other column is picked, and the pie, which cannot draw without one, says so. No scenario: the pie has no such message, and `PieChartLook.auto`, run on every removal, picks another categorical column. Needs the message in the pie (`validateState`) before a row can claim it. A layout whose category column the table lacks is the same case.

4. A Segment Angle Column holding NaN or ±Infinity leaves that value out, as an empty one (viewers treat such values as empty): right as it is (2026-10-30), requested below as an API test.

**API or Dart tests to write**

1. Dart regTest: more than 200000 categories shows "too many categories" (`pie_chart_core.dart`, `validateState`).
2. ApiTests: on a one-row table (demog-1000, first row) a pie of RACE draws 1 slice holding 100% of the circle.
3. ApiTests: with Segment Angle Column = WEIGHT (sum), one WEIGHT value set to NaN, Infinity or −Infinity leaves 4 slices and no error message.
4. ApiTests: the category column renamed (RACE to ORIGIN) is followed by the Category setting, no error.
5. ApiTests: with every row filtered out the pie draws 0 slices and reports no error; with the filter reset it draws its slices again.
6. ApiTests: once every row of the table is deleted the pie draws 0 slices and reports no error.
7. ApiTests: Max Radius set to 600 through the API (`maxRadius`, not typed in the panel: that path is item 2 of Decided) keeps 150: a known failure, core#1200, still open.

Moved out of BDD (no user gesture); requested from ApiTests.
