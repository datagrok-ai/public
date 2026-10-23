# Correlation plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `correlation_plot/CLAUDE.md` "Automation
surface" and the TestTrack `correlation-plot.md` / `correlation-plot-ui.md`. What could be said is in
`correlation-plot-colors-and-chrome.feature`.

## Waiting for the product

- **Selection highlighting in the diagonal histograms** (`correlation-plot-ui.md`, "Selection
  highlighting"; anchor "Diagonal histograms repaint"). The diagonal cells are built with the
  selection (`correlation_plot_core.dart:55`, `ColumnHistogram(..., dataFrame.selection)`, drawn in
  `#FF8C00` at alpha 150), but nothing repaints the matrix when the selection changes: on demog-1000
  selecting 300 rows leaves the canvas byte-identical (18 orange-ish pixels before and after, a hover
  changes nothing), and the orange bars appear only after an unrelated repaint — a window resize
  (548 px) or a filter change. The scenario is ready and passes the repaint claim the day the plot
  listens to `onSelectionChanged`:

  ```gherkin
  Scenario: A selection is drawn into the diagonal histograms and leaves the coefficients alone
    Then the "cell AGE x AGE" area of correlation plot viewer should not contain the color "#FF8C00"
    When user selects the first 20 rows
    Then 20 rows should be selected
    And the "cell AGE x AGE" area of correlation plot viewer should have repainted
    And the "cell AGE x AGE" area of correlation plot viewer should contain the color "#FF8C00"
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    When user clears the row selection
    Then the "cell AGE x AGE" area of correlation plot viewer should have repainted
    And the "cell AGE x AGE" area of correlation plot viewer should not contain the color "#FF8C00"
    And no errors should have been logged
  ```

## Needs a reading

- **Grid > Grid Color Coding > Apply to text** (GROK-19052). The coloured text is drawn on the inner
  grid's canvas; neither the plot's status nor the re-exported grid readings say whether a cell's
  text is coloured, so the scenario would end on a lone `no errors should have been logged`. Needs a
  `text color of cell <X> x <Y>` reading (the colour the cell's text was drawn with, from the grid's
  `GridCellStyle.textColor` at render time, `correlation_status.dart` beside `color of cell`):

  ```gherkin
  Scenario: Apply to text paints the coefficients instead of the cells
    When user picks "Grid > Grid Color Coding > Apply to text" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the "text color of cell HEIGHT x AGE" reading of correlation plot viewer should be "#c3c3ff"
    And the "color of cell HEIGHT x AGE" reading of correlation plot viewer should be "#ffffff"
    When user picks "Grid > Grid Color Coding > Apply to text" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the "color of cell HEIGHT x AGE" reading of correlation plot viewer should be "#c3c3ff"
    And no errors should have been logged
  ```

- **No cell tooltip while Tooltip > Visible is off** (`correlation-plot-ui.md`, "Tooltip"). The
  cell tooltip comes from the inner grid's hover, debounced by 200 ms
  (`grid/features/grid_tooltip.dart:6`), so `tooltip should not contain the text "Pearson R"` right
  after the hover passes before the tooltip could have appeared — it would pass with the menu item
  doing nothing. `correlation-plot-colors-and-chrome.feature` keeps the positive half (Visible off
  flips `showTooltip`, Visible on brings the tooltip back). Needs a tooltip-settled signal: a
  `tooltip pending` reading (or a render-pending flag the settle step waits on) set from the cell
  mouse-enter until the debounced handler has run and decided, in `grid_tooltip.dart`, so a
  negative claim is read only after the grid had its chance to show the tooltip:

  ```gherkin
  Scenario: With Tooltip > Visible off, a hovered cell shows no tooltip
    When user picks "Tooltip > Visible" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then "showTooltip" property of correlation plot viewer should be "false"
    When user clicks on the "cell STARTED x AGE" area of correlation plot viewer
    And user hovers over the "cell WEIGHT x HEIGHT" area of correlation plot viewer
    Then the grid tooltip should have settled
    And tooltip should not contain the text "Pearson R"
    When user moves the pointer away from correlation plot viewer
    And user picks "Tooltip > Visible" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then "showTooltip" property of correlation plot viewer should be "true"
    And no errors should have been logged
  ```
