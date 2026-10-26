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

## Edge cases (`radar-edge-cases.feature`)

Written for `radar-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Steps requested**

1. The Charts package on the local stand reports no Radar readings (`axes`, `rows shown`, `message`), so the feature claims the Values property the removal rewrites and the message text. Once the package with `radar-status.ts` is on the stand, the removal scenario should claim `the "axes" reading of radar viewer should be "AGE, WEIGHT"` and the empty row sets `the "rows shown" reading … should be 0`.

**Questions for owners**

1. The Min/Max values map with min above max for an axis: every value is clamped to max and the axis goes flat. Refuse or swap?
2. A Values column holding NaN or ±Infinity under Column normalization: how are the non-finite values drawn?
3. Show Current Row on a table with no current row highlights row 1 (`Math.max(currentRowIdx, 0)`, already in this file): should it highlight nothing?

**API or Dart tests to write**

1. Charts package test: a constant Values column under Column normalization gives ECharts a zero range.

### Won't do
