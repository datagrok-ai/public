# Line chart — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `line_chart/CLAUDE.md` "Automation
surface". What could be said about the mouse-over row line and the overview strip is in
`line-chart-mouse-over-and-overview.feature`.

## Needs a reading

- **Show Mouse-over Category** (`line-chart.md`, the Mouse-over group of the context menu). The
  property works: with the chart split by the column of a category hovered in another viewer (a bar
  of a bar chart on the same column), the line of that category is darkened and the others drawn at
  a fifth of their alpha while the flag is on, and nothing changes while it is off (walked by hand,
  23 Oct 2026; `showCat`, `line_chart_core.dart:97`, `renderers/line_chart_line_renderer.dart:245-251`).
  Nothing in `line_chart_status.dart` says which line is highlighted, so a claim could only be an
  ink diff — which would also move for the mouse-over row group.

  Needs a `mouse-over category` reading (the hovered category's split value, `''` when none or when
  `showCat` is false) beside `categories` in `line_chart_status.dart`:

  ```gherkin
  Scenario: A category hovered in a bar chart highlights its line, and Show Mouse-over Category turns that off
    Given user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName       | AGE    |
      | yColumnNames      | HEIGHT |
      | splitColumnNames  | SEX    |
    And user adds a bar chart viewer with:
      | splitColumnName | SEX |
    When user hovers over the "bar F" area of bar chart viewer
    Then the "mouse-over category" reading of line chart viewer should be "F"
    When user moves the pointer away from bar chart viewer
    And user sets "showMouseOverCategory" property of line chart viewer to "false"
    And user hovers over the "bar F" area of bar chart viewer
    Then the "mouse-over category" reading of line chart viewer should be ""
    When user moves the pointer away from bar chart viewer
    Then no errors should have been logged
  ```

## Edge cases (`line-chart-edge-cases.feature`)

Written for `line-chart-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Steps requested**

1. Regression line and moving average after a table swap (core#1100, wrong behaviour 1): `regression lines` and `moving average lines` read what the status reports, which may survive the dead subscriptions. Needs a reading of what the last frame actually drew (the draw calls of the regression and moving-average passes), then: `When user sets "table" property … to "demog-1000"` and back, `Then the "regression lines drawn" reading of line chart viewer should be 1`, tagged `@known-failure … @realizes:core#1100`.
2. Stacked Area under a filter (core#1100, wrong behaviour 6): the expected Y maximum is the stack of the rows shown. Needs a fixture with a known stacked maximum per filter state (a small `user opens a table … with:` table whose stacked maximum is computed by hand), then `the 'y axis max of "<col>"' reading … should be between …` tagged with core#1100.

**Questions for owners**

1. Lower Control Limit above Upper Control Limit (e.g. 5 and 0): swap, refuse, or report? `line_chart_spc_feature.dart` has no cross-check.
2. X Min above X Max (and Y Min above Y Max): the line chart swaps them silently, the histogram reports "X axis Max should be greater than Min". Which is the rule?
3. X Min / X Max outside every X value (700000..800000 on spgi-100): the setting is kept and the view does not move. Intended?
4. Line Width 0 or negative: the canvas ignores it and keeps the previous width. Should the property have a minimum?
5. The chart bound to another table whose view is then closed: rebind to the view's table, or show an empty chart with a message?
6. A layout whose X, Y and split columns the table does not have: `auto()` re-picks X and Y and drops the split silently. Re-pick or report the missing columns?

**API or Dart tests to write**

1. Dart regTest: a layout or `setOptions` with `chartTypes: ['Bar Chart']` or an unknown type throws in `LineChartRenderer.getRenderer` on every frame.
2. Dart regTest: the viewport-restore timer armed by `afterDeserialization` fires on a chart closed within 4 s (core#1100).
3. Dart regTest: the aggregated marker map selects (61, Other) for the rows of (61, R) (core#1100, wrong behaviour 7); hidden rows hit-tested under Row Source Selected (wrong behaviour 2).

### Won't do

- Whiskers over single-row groups, ±Inf in X and Y, all-null columns: core#1100 found them clean and the crash-only Data Manipulation test covers the class.
- The SPC "Invalid formula" path of the viewer Filter: written once for three viewers in `row-source/row-source-edge-cases.feature`.
