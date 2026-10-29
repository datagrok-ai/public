# Radar edge cases — what the features cannot say yet

Written for `radar-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Steps requested**

1. The Charts package on the local stand reports no Radar readings (`axes`, `rows shown`, `message`), so the feature claims the Values property the removal rewrites and the message text. Once the package with `radar-status.ts` is on the stand, the removal scenario should claim `the "axes" reading of radar viewer should be "AGE, WEIGHT"` and the empty row sets `the "rows shown" reading … should be 0`.

**Decided, steps requested** (the rules of `viewer-rules-edge-cases.feature`)

1. The Min/Max values map with min above max for an axis is swapped. Today every value is clamped to max and the axis goes flat. The radar reports no reading of an axis's drawn range, so no row claims it yet.
2. Show Current Row on a table with no current row highlights nothing. Today it highlights row 1 (`Math.max(currentRowIdx, 0)`). The radar reports no reading of the highlighted row.
3. A row with NaN or ±Infinity in a Values column under Column normalization draws no point for that value. The feature claims the rest (the chart painted, the Values property, no error logged); it does not claim the missing point, because the radar reports no reading of the points a row line draws. Requested: a reading of the drawn points per row line (e.g. `row N points`: the axes the line has a point on). NaN reaches the table only as an empty cell (a CSV `NaN`, and a formula's `Sqrt(-1)`, both read as missing), so the NaN half is the missing-value case.

**API or Dart tests to write**

1. Charts package test: a constant Values column under Column normalization gives ECharts a zero range.

## Won't do

Nothing.
