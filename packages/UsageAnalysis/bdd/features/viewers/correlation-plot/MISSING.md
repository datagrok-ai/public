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

## Edge cases

Written in the edge-case round of 2026-10-26 (localhost). No library, binding or code change was made. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Decided** (2026-10-28): a column removed from narrowed axes leaves only that column: core#1172,
   fixed. No row passing, Row Source Selected with one row and a table with no numeric column give an
   empty plot with no message; that is right as it is.
2. **API or Dart tests to write**: a table with 120 numeric columns opens with Show Coefficients off.
   On demog-1000 with the plot open:
   - AGE renamed to AGE2: both axes read "AGE2, HEIGHT, WEIGHT, STARTED", the cell HEIGHT x AGE2
     exists (HEIGHT x AGE does not) and shows -0.23, with no error.
   - Axes narrowed to AGE, HEIGHT, WEIGHT, then WEIGHT removed: both axes read "AGE, HEIGHT", with
     no error (core#1172).

   Moved out of BDD (no user gesture); requested from ApiTests.
