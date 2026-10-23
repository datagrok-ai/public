# Box plot — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the box plot's "Automation surface" (`core/client/d4/lib/src/viewers/box_plot/CLAUDE.md`)
and `files/TestTrack/Viewers/BoxPlot/box-plot-ui.md`. What could be said went into
`box-plot-hover-and-date-categories.feature` (the row-group hover and Show Mouse Over Row Group,
github-3065, the month and quarter categories) and into claims added to `box-plot.feature`
(category selector), `box-plot-selection.feature` (Show Mouse Over Point's positive half),
`box-plot-group-comparison.feature` (the Adjust by selector goes) and `box-plot-filter.feature`
(the colour scale keeps clear of both axes, GROK-20171).

## Still open

- **The mouse-over highlight follows the pointer** (`box-plot-ui.md` "Mouse-over highlight", steps
  2-4). The status reports one `marker` area — a single drawn marker — and nothing about which row
  the pointer is over, so "move to a different point — the highlight follows" has no second point to
  aim at and no reading to compare. Needs either a second marker area (`marker of row <r>` for the
  markers drawn, as the scatter plot reports `marker of row <r>`) or a `hovered row` reading (from 1,
  0 for none, as the scatter plot publishes), from the box plot's mouse-over hit test in
  `box_plot_core.dart`. Then:

  ```gherkin
  Scenario: The mouse-over highlight follows the pointer from marker to marker
    When user hovers over the "marker of row 11" area of box plot viewer
    Then the "hovered row" reading of box plot viewer should be 11
    And the mouse-over row of the table should be 11
    When user hovers over the "marker of row 12" area of box plot viewer
    Then the "hovered row" reading of box plot viewer should be 12
    When user moves the pointer away from grid
    Then the "hovered row" reading of box plot viewer should be 0
    And no errors should have been logged
  ```
