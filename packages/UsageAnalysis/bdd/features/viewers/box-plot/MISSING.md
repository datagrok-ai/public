# Box plot — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the box plot's "Automation surface" (`core/client/d4/lib/src/viewers/box_plot/CLAUDE.md`)
and `files/TestTrack/Viewers/BoxPlot/box-plot-ui.md`. What could be said went into
`box-plot-category-label-hover.feature` (the row-group hover and Show Mouse Over Row Group) and
into claims added to `box-plot.feature` (category selector, the month and quarter categories),
`box-plot-selection.feature` (Show Mouse Over Point's positive half, github-3065),
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

- **The other groups fade under a cross-viewer hover** (`boxplot-pointer-select.md` Scenario 7, step 4,
  github-3065). `box-plot-selection.feature` claims the box plot repaints and the
  hovered Caucasian group keeps its own RACE colour, not that Asian, Black and Other fade. The fade is
  a translucent back-colour overlay over the whole plot (`fillRect(ctx, overlay.client,
  alpha(look.backColor, 140))` in `box_plot_core.dart`), and the colour step's tolerance still finds
  Asian's #1F77B4 under it (a run of 23 Oct: `the "Asian values" area ... should not contain the color
  "#1F77B4"` failed with 6469 px matched while the screenshot shows the group faded). Needs a
  `faded groups` reading (the categories drawn under the overlay in the last pass, `""` when none) or
  an `<category> faded` flag next to `<category> values` in the box plot's status. Then:

  ```gherkin
  Scenario: A hover on another viewer's bar fades the other groups
    When user adds a bar chart viewer with:
      | Split | RACE |
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then the "faded groups" reading of box plot viewer should be "Asian, Black, Other"
    When user moves the pointer away from grid
    Then the "faded groups" reading of box plot viewer should be ""
    And no errors should have been logged
  ```

## Edge cases (`box-plot-edge-cases.feature`)

Added with `box-plot-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Steps requested**:
   - The log-axis warning: the box plot draws "… can't be visualized in the log scale" over a value
     with zero or negative values but reports no area for it. Wanted: a `y warning` area as the density
     plot reports.
   - Alpha typed in the property panel as 0, 1, 1.5 or -0.1 (help: 0 < alpha < 1) needs a step that
     types into a property row (the property step goes through the API).
2. **Data**: the group comparison shown at 30 categories and hidden at 31 (help box-plot.md) needs a
   categorical column with exactly 30 and 31 categories; demog-1000 has none.
3. **Questions for owners**: the Category 1 column removed while in use; a group of one row (does the
   info tooltip give the variance reason); the control group filtered out.
4. **Dart tests to write**: a value column holding NaN (not an empty value) with the violin style, and
   retained listeners after ten add/close cycles and ten table swaps (core#960).
