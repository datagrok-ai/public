# Pie chart edge cases — what the features cannot say yet

Written for `pie-chart-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Decided** (`viewer-rules-edge-cases.feature`)

1. A renamed category column is followed (a row of the rename outline); every row filtered out is a row of the no-rows outline.
2. Max Radius typed above 500 is not taken (a known failure in the out-of-bounds outline); Start Angle typed as −30 or 400 is the same, written in `pie-chart-edge-cases.feature`. Both are core#1200.
3. The category column removed: no other column is picked, and the pie, which cannot draw without one, says so. No scenario: the pie has no such message, and `PieChartLook.auto`, run on every removal, picks another categorical column. Needs the message in the pie (`validateState`) before a row can claim it. A layout whose category column the table lacks is the same case.

4. A Segment Angle Column holding NaN or ±Infinity leaves that value out, as an empty one (viewers treat such values as empty): right as it is (2026-10-30), claimed in `pie-chart-edge-cases.feature` (sum of WEIGHT with one value replaced; 4 slices, no message).

**API or Dart tests to write**

1. Dart regTest: more than 200000 categories shows "too many categories" (`pie_chart_core.dart`, `validateState`).
