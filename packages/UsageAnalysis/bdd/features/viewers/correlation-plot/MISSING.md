# Correlation plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `correlation_plot/CLAUDE.md` "Automation
surface" and the TestTrack `correlation-plot.md` / `correlation-plot-ui.md`. What could be said is in
`correlation-plot-colors-and-chrome.feature`.

## Waiting for the product

- **Selection highlighting in the diagonal histograms** is now a scenario in
  `correlation-plot-colors-and-chrome.feature` (core#1087, fixed).

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

## Edge cases (`correlation-plot-edge-cases.feature`)

Added with `correlation-plot-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Questions for owners**: axes narrowed to AGE, HEIGHT, WEIGHT and then WEIGHT removed (today both
   axes reset to every numeric column); what the cells show with no row passing; Row Source Selected
   with a single row; a table with no numeric column.
2. **ApiTests to write**: a table with 120 numeric columns opens with Show Coefficients off.
