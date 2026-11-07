# Radar — what the features cannot say yet

Written during the Charts gap round of 23 Oct 2026, against the viewer's "Automation surface"
(`public/packages/Charts/src/viewers/radar/CLAUDE.md`), the old spec `Charts/playwright/radar.test.ts`
(step 10) and the library's rules (`public/libraries/bdd/CLAUDE.md`). `radar-current-row.feature`
holds what could be said: Show Current Row on draws a `#00FF00` line inside the `view` area, off
leaves none.

## Still open

- **Which row the highlighted line belongs to.** The colour claim proves a current-row line is drawn,
  not that it is the current row's: a highlight stuck on row 1 (`Math.max(currentRowIdx, 0)` in
  `init()`) would satisfy it after any grid click. Needs a `current row` reading in
  `Charts/src/viewers/radar/radar-status.ts` — the 1-based row whose line the frame drew with the
  current-row colour (series 2, the item with a symbol and `itemStyle.color === currentRowColor`), or
  `0` when Show Current Row is off. Then:

  ```gherkin
  Scenario: The current row's line follows the grid's current row
    When user clicks on the "cell 5 of AGE" area of grid
    Then row 5 should be current
    And the "current row" reading of radar viewer should be 5
    When user clicks on the "cell 9 of AGE" area of grid
    Then the "current row" reading of radar viewer should be 9
    When user clicks on settings icon of radar viewer
    Given "Selection" category in context panel is expanded
    When user unchecks "Show Current Row" property in context panel
    Then the "current row" reading of radar viewer should be 0
    And no errors should have been logged
  ```

- **Show Mouse-Over Row.** The mouse-over line is drawn in `#AAAAAA`, a grey the axis and split-line
  ink also contain, so a colour claim cannot tell it apart. Needs a `mouse-over row` reading from the
  same series (the row drawn with `mouseOverRowColor`), and the scenario hovers `cell 7 of AGE`, claims
  `mouse-over row` 7, unchecks Show Mouse-Over Row, hovers again and claims 0.
- **A line's hit area.** The chart reports no area per line (`line <row>`), so clicking a line to make
  its row current (the viewer's `click` handler sets `currentRowIdx`) has no gesture.

## Edge cases

The edge cases of the Radar viewer that no feature states, with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Steps requested**

1. The Charts package on the local stand reports no Radar readings (`axes`, `rows shown`, `message`). Once the package with `radar-status.ts` is on the stand, a check can claim the axes (`the "axes" reading of radar viewer should be "AGE, WEIGHT"`) and an empty row set (`the "rows shown" reading … should be 0`).

**Decided, steps requested** (the rules every viewer shares, listed in `viewers/MISSING.md`)

1. The Min/Max values map with min above max for an axis is swapped. Today every value is clamped to max and the axis goes flat. The radar reports no reading of an axis's drawn range, so no row claims it yet.
2. Show Current Row on a table with no current row highlights nothing. Today it highlights row 1 (`Math.max(currentRowIdx, 0)`). The radar reports no reading of the highlighted row.
3. A row with NaN or ±Infinity in a Values column under Column normalization draws no point for that value; the rest of the chart draws. No check claims the missing point, because the radar reports no reading of the points a row line draws. Requested: a reading of the drawn points per row line (e.g. `row N points`: the axes the line has a point on). NaN reaches the table only as an empty cell (a CSV `NaN`, and a formula's `Sqrt(-1)`, both read as missing), so the NaN half is the missing-value case.

**API or Dart tests to write**

1. Charts package test: a constant Values column under Column normalization gives ECharts a zero range.
2. HEIGHT removed from the table: Values is `AGE, WEIGHT`, the chart is painted, no column message, no error; with AGE and WEIGHT removed too, Values is empty and "The Radar viewer requires a minimum of 1 numerical column." is shown.
3. Every row filtered out (AGE between 200 and 300): the radar redraws, shows no column message and logs no error.
4. Every row of the table deleted: the radar redraws, shows no column message and logs no error.
5. Values columns holding NaN, Infinity and -Infinity under Column normalization: Values is `a, b, c`, no column message, the chart is painted, no error.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

Nothing.
