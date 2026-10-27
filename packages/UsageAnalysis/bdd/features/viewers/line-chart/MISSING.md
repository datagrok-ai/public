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

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed line chart issues into `line-chart-row-source-regressions.feature`
(the split legend ones are `line-chart-regressions.feature`). Each item needs a step or a reading
that does not exist yet.

- **The "too many split combinations" message (public#3288) is gone — nothing to claim** (Olesia,
  27 Oct 2026). 7c24a9df21 (18 Mar 2025) warned "Selected split columns resulted in N combinations.
  Supported limit is …"; "Line chart: Optimizations and speedups" (b374989c2c, the next day) commented
  it out and the d4 cleanup c0d2266d64 (Sep 2026) deleted it. Its limit was 100,000,000 combinations
  (`TOTAL_CATEGORIES_LIMIT_TO_STORE`), out of reach on demog anyway. What public#3210 asked for — a
  split this large does not break the page — is claimed in `line-chart-row-source-regressions.feature`.

- **A click on a line selects its category; Control adds, Control+Shift removes (public#2608,
  public#2360).** The chart hit-tests a split line (`LineChartCore.hitTestLineCategory`,
  `line_chart_core.dart:808`, through `renderers/line_chart_line_renderer.dart:497`) and a click on it
  selects the category's rows — under Row Source = Selected it used to clear the selection instead
  (public#2608). `line_chart_status.dart` reports `point of row <n>` but no area on a line, so there is
  nothing to click that is known to be on a segment and not on a marker. Needs a `line <category>` area
  (a point on a drawn segment, away from the vertices, that `hitTestLineCategory` accepts — as the
  histogram reports `line <category> point`):

  ```gherkin
  Scenario: Under Row Source Selected a click on a line selects that line's rows (public#2608)
    Given user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName      | AGE      |
      | yColumnNames     | HEIGHT   |
      | splitColumnNames | RACE     |
      | rowSource        | Selected |
    When user selects rows where "RACE" is one of "Asian, Black"
    Then the "lines" reading of line chart viewer should be 2
    When user clicks on the "line Asian" area of line chart viewer
    Then only rows where "RACE" is "Asian" should be selected
    And the "lines" reading of line chart viewer should be 1
    And no errors should have been logged

  Scenario: A line click selects its category, Control adds one, Control+Shift takes it away (public#2360)
    Given user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName      | AGE    |
      | yColumnNames     | HEIGHT |
      | splitColumnNames | RACE   |
    When user clicks on the "line Asian" area of line chart viewer
    Then only rows where "RACE" is "Asian" should be selected
    When user clicks on the "line Black" area of line chart viewer holding Control
    Then only rows where "RACE" is one of "Asian, Black" should be selected
    When user clicks on the "line Asian" area of line chart viewer holding Control+Shift
    Then only rows where "RACE" is "Black" should be selected
    And no errors should have been logged
  ```

- **A split or a Y column removed on the plot shows up in the property panel at once (public#2951,
  public#2875, public#2623, public#2880).** The on-plot Y and split selectors carry a cross per
  column; no step clicks the cross of one named column (`user adds {string} to the splits of {widget}`
  exists, its removal does not), and the property panel's `Split` / `Y` rows are read only through
  their text. Needs `user removes the {string} split of {widget}` and `user removes the {string} Y
  column of {widget}` (the cross of that column's chip in the inner box selector), then:

  ```gherkin
  Scenario: A split removed on the plot leaves the property panel at once (public#2951)
    Given user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName  | AGE            |
      | yColumnNames | HEIGHT, WEIGHT |
    When user clicks on settings icon of line chart viewer
    And user adds "SEX" to the splits of line chart viewer
    Then the "split columns" reading of line chart viewer should be 1
    And "Split" property in context panel should contain text "SEX"
    When user removes the "SEX" split of line chart viewer
    Then the "split columns" reading of line chart viewer should be 0
    And "Split" property in context panel should not contain text "SEX"
    When user removes the "WEIGHT" Y column of line chart viewer
    Then the "y columns" reading of line chart viewer should be "HEIGHT"
    And "Y" property in context panel should not contain text "WEIGHT"
    And no errors should have been logged
  ```

- **A logarithmic axis over non-positive values shows its warning (public#3701).** The chart draws
  a warning next to such an axis (`_renderLogWarnings`, `line_chart_core.dart:2176`, boxes
  `_xLogWarningBox` and the per-column Y boxes), but the status reports neither box. The density plot
  reports the same warnings as `x warning` / `y warning` (`density_plot_status.dart:8-9`). Needs the
  same two areas in `line_chart_status.dart`, present only while the warning is drawn:

  ```gherkin
  Scenario: A logarithmic Y axis over values below zero shows its warning, and a linear one drops it (public#3701)
    Given user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | Average Mass     |
      | yColumnNames | Chemical Space X |
    When user sets "yAxisType" property of line chart viewer to "logarithmic"
    Then line chart viewer should have a "y warning" area
    When user sets "yAxisType" property of line chart viewer to "linear"
    Then line chart viewer should not have a "y warning" area
    And no errors should have been logged
  ```

- **After a category filter the X labels stay over their points (public#1456).** A categorical X axis
  reports no label areas, so "the label of Chemist X sits over the point of its row" has nothing to
  compare. Needs `x label <category>` areas (the drawn label's box) on the categorical X axis; the
  claim is then horizontal alignment with `point of row <n>`, which needs a geometry step on two
  areas (`the {string} area of {widget} should be centred over the {string} area`):

  ```gherkin
  Scenario: Filtered categories keep their labels over their points (public#1456)
    Given user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | Chemist 521      |
      | yColumnNames | Chemical Space Y |
    When user filters rows where "Chemist 521" is one of "<a>, <c>, <e>"
    Then the "x label <c>" area of line chart viewer should be centred over the "point of row <r>" area
    And no errors should have been logged
  ```
  (`<a>`, `<c>`, `<e>` and the row `<r>` of a `<c>` point to be read from the data when the step exists.)

- **A custom category order is what the axis and the points follow (public#1551).** No step sets a
  column's custom category order (the `.category-order` tag, set through Sort > Custom in the column
  menu), and the line chart reports no area for a point's X label, so "the tooltip's X is the label
  under it" cannot be put. Needs `user orders the categories of {string} column as {string}` and the
  `x label <category>` areas above:

  ```gherkin
  Scenario: Under a custom category order a point's tooltip names the category under it (public#1551)
    Given user opens spgi dataset
    When user orders the categories of "Stereo Category" column as "S_UNKN, R_ONE, S_ABS, S_ACHIR, S_PART"
    And user adds a line chart viewer with:
      | xColumnName  | Stereo Category |
      | yColumnNames | Average Mass    |
    Then the "x label S_UNKN" area of line chart viewer should lie to the left of the "x label R_ONE" area
    When user hovers over the "point of row 1" area of line chart viewer
    Then the tooltip should show "Stereo Category" as "<category of row 1>"
    And no errors should have been logged
  ```

- **Multi Axis with Y Global Scale puts a constant series on its own tick (core#823).** The fix is a
  matter of a few pixels; `y axis min/max of "<col>"` give the scale, not where a line was drawn.
  Needs a reading of the value a drawn point maps back to (`y of point of row <n> of "<col>"`, from
  the point's screen Y through the chart's transform), then:

  ```gherkin
  Scenario: Under Multi Axis and Y Global Scale a constant series sits on its tick (core#823)
    Given user opens a table "flat" with:
      | x | a  | b  | c  |
      | 0 | 10 | 20 | 30 |
      | 1 | 10 | 20 | 30 |
      | 2 | 10 | 20 | 40 |
    And user adds a line chart viewer with:
      | xColumnName  | x       |
      | yColumnNames | a, b, c |
      | multiAxis    | true    |
      | yGlobalScale | true    |
    Then the "y of point of row 1 of \"a\"" reading of line chart viewer should be between 9.95 and 10.05
    And the "y of point of row 1 of \"b\"" reading of line chart viewer should be between 19.95 and 20.05
    And no errors should have been logged
  ```

- **Resizing does not join the first and last points (public#1671), and the SPC average label keeps
  off the axis labels (public#4013).** Neither has a reading: the first needs the number of segments
  a series drew (`segments of "<col>"`, unchanged by a resize), the second an `average label` area
  next to the `average` band, to be claimed with `should not overlap` against `y axis`.

- **A split column with too many categories (public#3288, public#3210).** The only limit found in
  the code is the split picker's preview gate (`maxCategoriesPreview: 200` in `line_chart.g.dart`,
  `column_combo_box.dart:509`); the "not added: too many categories" message the issue's screenshot
  shows was not found. To be walked by hand first — what a person sees after picking `Id` (100
  categories in spgi-100; the full spgi has more) — and claimed by its balloon or text once known.
