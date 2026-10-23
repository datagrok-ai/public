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
