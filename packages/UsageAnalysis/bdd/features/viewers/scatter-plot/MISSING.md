# Scatter plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `scatterplot/CLAUDE.md` "Automation
surface", the TestTrack `scatterplot-*.md` and `scatter-plot-ui.md`. The "(no value)" legend click is
in `scatter-plot-legend-no-value.feature`; the rejected line formula (GROK-17282) in
PowerPack's `formula-lines-dialog.feature`.

## Needs a reading

- **Structures on the X axis** (`scatter-plot-ui.md`, "Structures as markers and as an axis", items
  7-8, spgi-100). With a molecule column on X the axis draws structures instead of tick labels;
  the status reports the axis column and range, not what the ticks are drawn as. Needs an
  `x axis labels drawn as` reading (`text` / `structure`, from the axis painter's renderer), claimed
  after the pick and again after a saved layout is applied.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed scatter plot issues into `scatter-plot-regressions.feature`,
`scatter-plot-look-regressions.feature` and `scatter-plot-labels-regressions.feature`.

- **The missing-column message (public#3091).** A layout applied after its Y column was removed draws
  `y: Column "HEIGHT" does not exist.` on the plot (`DataFrameViewer._validateColumns`), but
  `_ScatterPlotStatus.of` builds its `WidgetStatus` without `error:` (the line chart passes
  `error: plot.error`), so `scatter plot viewer should report the error …` reads `""`. Removing the
  column live is no test either: the plot re-picks another column. Needs `error: p.error` in
  `scatterplot_status.dart`:

  ```gherkin
  Scenario: A plot whose Y column is gone says which column it misses (public#3091)
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should report no error
    When user saves the layout of the current table view
    And user removes "HEIGHT" column
    And user loads the saved layout
    Then scatter plot viewer should report the error "y: Column \"HEIGHT\" does not exist."
    And no errors should have been logged
  ```
- **Markers the renderer drew (public#2108, public#1411).** `rows shown` counts the rows the axes can
  place (`isRowDrawable`), not what the renderer drew, so it reads 872 whether the dot renderer puts the
  128 rows with no HEIGHT on the zero line (public#2108) and whether a numeric Color column's markers
  come back after a layout (public#1411, which said markers went missing). A pixel claim does not
  help: a layout round trip repaints the canvas by a few pixels with every marker in place (measured
  26 Oct 2026). Needs a `markers drawn` reading counted by `_ScatterPlotRenderer` in the last frame:

  ```gherkin
  Scenario: Dot markers draw no row whose Y is empty (public#2108)
    When user sets "Marker Type" property of scatter plot viewer to "dot"
    Then the "markers drawn" reading of scatter plot viewer should be 872
    When user sets "Marker Type" property of scatter plot viewer to "circle"
    Then the "markers drawn" reading of scatter plot viewer should be 872
    And no errors should have been logged

  Scenario: A numeric Color column comes back from a layout with every marker it drew (public#1411)
    When user sets "Color" property of scatter plot viewer to "AGE"
    And user saves the layout of the current table view
    And user loads the saved layout
    Then the "markers drawn" reading of scatter plot viewer should be 872
    When user sets "Color" property of scatter plot viewer to ""
    Then no errors should have been logged
  ```
- **A wheel zoom ends after the plot says it has settled.** After `user scrolls the mouse wheel up 3
  times over the "view" area`, `user remembers the value range` and an immediate `should show the
  remembered value range` failed: the window went on narrowing (x 79.3..127.3, then 85.5..121.1), so
  the zoom animation is not in `isRenderPending`. The regression features set a window with the range
  slider instead; a claim on a wheel zoom's end needs the animation in the pending flags.
