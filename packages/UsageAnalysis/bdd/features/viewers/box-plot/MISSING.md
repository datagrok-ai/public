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

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed box plot issues into `box-plot-regressions.feature`.

- **The violin style over a value column holding a real NaN (core#960).** The crash needed NaN, not a
  missing value: a column built with `DG.Column.fromList`, where NaN is not none, threw `NaN.floor()` in
  `_defineCatKdeValues` and on every hover. Every data step writes NaN as a missing value: `user sets
  "WEIGHT" column in row 2 to "NaN"` stores the float null (`isNone` true, measured 26 Oct 2026), and so
  do the calculated columns `Sqrt(-1)`, `0/0` and `${WEIGHT} * NaN`. Needs a step that writes a raw NaN
  into a float column's buffer (`col.getRawData()[i] = NaN` followed by `fireValuesChanged`), then:

  ```gherkin
  Scenario: The violin style draws and hovers over a value column holding a NaN (core#960)
    Given user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | WEIGHT |
      | Category 1 | RACE   |
      | Plot Style | violin |
    When user writes a raw NaN into "WEIGHT" column in row 2
    And user hovers over the "Caucasian values" area of box plot viewer
    Then box plot viewer should report no error
    And box plot viewer should have a "category Caucasian" area
    And no errors should have been logged
  ```

- **The statistics strip lines up under each box and uses the column's format (public#3693,
  public#3773).** The status reports one `stats` area for the whole strip and no per-category cell or
  text. Needs `stats <category>` areas (one per category column of the strip) and a `stats <category>`
  reading of the text drawn there:

  ```gherkin
  Scenario: Each category's statistics sit under its own box in the column's format (public#3693)
    Given user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Category 1      | RACE |
      | Value           | AGE  |
      | Show Statistics | true |
    Then the "stats Caucasian" area of box plot viewer should lie below the "Caucasian values" area
    And the "stats Caucasian" area of box plot viewer should lie to the right of the "stats Black" area
    And no errors should have been logged
  ```
