# Line chart — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `line_chart/CLAUDE.md` "Automation
surface". What could be said about the mouse-over row line and the overview strip is in
`line-chart-mouse-over-and-overview.feature`.

## Needs a reading — and a question for the viewer's owner

- **Show Mouse-over Category** (`line-chart.md`, the Mouse-over group of the context menu). The
  property reaches the picture by two paths, and neither is reported:
  - The line renderer: `showCat` (`line_chart_core.dart:97`) is true when the chart is split by the
    column of the category hovered in another viewer (`AppEvents.mouseOverCategoryColumn`, e.g. a
    bar of a bar chart on the same column) and the property is on; the line of that category is
    then darkened and the others drawn at a fifth of their alpha
    (`renderers/line_chart_line_renderer.dart:245-251`). Nothing in `line_chart_status.dart` says
    which line is highlighted, so a claim could only be an ink diff — which would also move for the
    mouse-over row group.
  - The overlay: `if (look.showMouseOverCategory && _catDraw)` (`line_chart_core.dart:1977`) dims
    the chart and redraws it, but the static `_catDraw` is only ever reset to `false` (line 1984) and
    never set to `true`, so this branch never runs. **Candidate product question**: is the overlay
    dead code, or was a `_catDraw = true` lost when the category hover moved to `showCat`? Not
    filed — for the owner to answer.

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
