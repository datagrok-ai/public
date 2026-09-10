@journey @viewers @realizes:viewers.correlation-plot
Feature: Correlation plot — the matrix and the numbers in it
  What the matrix holds: one cell per pair of numerical columns, the coefficient behind each cell,
  the two decimals it draws there, the histogram on the diagonal, and what Show Pearson R and the
  correlation type change about all of that.
  The spec this replaces read its numbers by calling `cp.getCorrelation(c1, c2)` on the viewer
  object — nine times across three files — and never looked at the matrix at all; where it did look,
  it sampled five canvas pixels per cell with a median filter after a six-attempt probe-click loop
  that mutated its own `pinnedW`/`headerH` constants until a click landed. Both are gone: a cell is
  the region `cell HEIGHT x AGE`, its number is `correlation of HEIGHT and AGE`, and its drawn text
  is `text of cell HEIGHT x AGE`. Those are three different facts and the scenarios below keep them
  apart — Show Pearson R empties the text and leaves the coefficient exactly where it was.
  The coefficients are checked twice over: against the values measured on this fixture, and against
  `DG.Stats` computed independently over the rows the plot's Row Source names — the only check here
  that does not ask the viewer the same question twice.
  demog-1000 has 1000 rows of which 128 have a blank HEIGHT, and both `plot.correlation` and
  `DG.Stats` are pairwise-complete: AGE/HEIGHT is computed over 872 rows while AGE/WEIGHT uses all
  1000, so `rows shown` is 1000 for the whole matrix and neither number is over "the rows shown".

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then 1000 rows should pass the filter
    And the "rows shown" reading of correlation plot viewer should be 1000
    And the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And the "y columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And the "cells" reading of correlation plot viewer should be 16
    And the "correlation type" reading of correlation plot viewer should be "Pearson"
    And the "show pearson r" reading of correlation plot viewer should be "true"
    And correlation plot viewer should be painted

  Scenario: Both axes are the table's numerical columns, and only those
    Then the "numerical columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And the "x columns" and "numerical columns" readings of correlation plot viewer should be the same
    And the "y columns" and "numerical columns" readings of correlation plot viewer should be the same
    And correlation plot viewer should have a "cell HEIGHT x AGE" area
    And correlation plot viewer should not have a "cell SEX x AGE" area
    And correlation plot viewer should not have a "cell RACE x AGE" area
    And the "columns shown" reading of correlation plot viewer should be 6
    And the "column order" reading of correlation plot viewer should be "__t, __name, AGE, HEIGHT, WEIGHT, STARTED"
    And correlation plot viewer should have a "row header HEIGHT" area
    And correlation plot viewer should have a "type header HEIGHT" area
    And no errors should have been logged

  Scenario: Every off-diagonal cell holds its pair's coefficient, and DG.Stats agrees
    Then the correlation of "HEIGHT" and "AGE" of correlation plot viewer should match the Pearson coefficient of the table
    And the correlation of "WEIGHT" and "AGE" of correlation plot viewer should match the Pearson coefficient of the table
    And the correlation of "WEIGHT" and "HEIGHT" of correlation plot viewer should match the Pearson coefficient of the table
    And the correlation of "STARTED" and "AGE" of correlation plot viewer should match the Pearson coefficient of the table
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "correlation of WEIGHT and AGE" reading of correlation plot viewer should be between 0.0647 and 0.0649
    And the "correlation of WEIGHT and HEIGHT" reading of correlation plot viewer should be between 0.4124 and 0.4125
    And the "correlation of STARTED and AGE" reading of correlation plot viewer should be between -0.0090 and -0.0089
    And the "correlation of AGE and HEIGHT" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And no errors should have been logged

  Scenario: The diagonal draws a histogram instead of a coefficient of one
    Then the "cell type of AGE x AGE" reading of correlation plot viewer should be "histogram"
    And the "cell type of HEIGHT x HEIGHT" reading of correlation plot viewer should be "histogram"
    And the "cell type of HEIGHT x AGE" reading of correlation plot viewer should be "correlation"
    And correlation plot viewer should have a "cell AGE x AGE" area
    And the "text of cell AGE x AGE" reading of correlation plot viewer should be ""
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And no errors should have been logged

  Scenario: What the cell drew is the coefficient rounded to two decimals
    Then the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "text of cell WEIGHT x HEIGHT" reading of correlation plot viewer should be "0.41"
    And the "correlation of WEIGHT and HEIGHT" reading of correlation plot viewer should be between 0.4124 and 0.4125
    And the "text of cell STARTED x AGE" reading of correlation plot viewer should be "-0.01"
    And no errors should have been logged

  Scenario: Show Pearson R empties the cells and halves their width, and moves no coefficient
    Then the "cell width" reading of correlation plot viewer should be 40
    And the "column width of HEIGHT" reading of correlation plot viewer should be 40
    When user sets "showPearsonR" property of correlation plot viewer to "false"
    Then the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be ""
    And the "text of cell WEIGHT x HEIGHT" reading of correlation plot viewer should be ""
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "cell width" reading of correlation plot viewer should be 20
    And the "column width of HEIGHT" reading of correlation plot viewer should be 20
    And correlation plot viewer should have repainted
    When user sets "showPearsonR" property of correlation plot viewer to "true"
    Then the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And the "cell width" reading of correlation plot viewer should be 40
    And no errors should have been logged

  Scenario: Spearman puts a different coefficient in the same cells
    When user sets "correlationType" property of correlation plot viewer to "Spearman"
    Then the "correlation type" reading of correlation plot viewer should be "Spearman"
    And the correlation of "HEIGHT" and "AGE" of correlation plot viewer should match the Spearman coefficient of the table
    And the correlation of "WEIGHT" and "HEIGHT" of correlation plot viewer should match the Spearman coefficient of the table
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2400 and -0.2399
    And the "correlation of WEIGHT and AGE" reading of correlation plot viewer should be between 0.0943 and 0.0944
    And the "correlation of WEIGHT and HEIGHT" reading of correlation plot viewer should be between 0.4453 and 0.4454
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.24"
    And correlation plot viewer should have repainted
    When user sets "correlationType" property of correlation plot viewer to "Pearson"
    Then the "correlation type" reading of correlation plot viewer should be "Pearson"
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And no errors should have been logged

  Scenario: Narrowing the axes re-tiles the matrix to their product
    When user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT |
      | yColumnNames | AGE, HEIGHT         |
    Then the "cells" reading of correlation plot viewer should be 6
    And the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT"
    And the "y columns" reading of correlation plot viewer should be "AGE, HEIGHT"
    And the "columns shown" reading of correlation plot viewer should be 5
    And correlation plot viewer should have a "cell WEIGHT x HEIGHT" area
    And correlation plot viewer should not have a "cell STARTED x AGE" area
    And correlation plot viewer should not have a "cell AGE x WEIGHT" area
    And the correlation of "WEIGHT" and "HEIGHT" of correlation plot viewer should match the Pearson coefficient of the table
    When user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of correlation plot viewer should be 16
    And correlation plot viewer should have a "cell STARTED x AGE" area
    And no errors should have been logged

  Scenario: A column with no variance has no coefficient and its cells stay blank
    Given user adds a calculated column "FLAT" with formula "0"
    When user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED, FLAT |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED, FLAT |
    Then the "cells" reading of correlation plot viewer should be 25
    And the "cell type of FLAT x AGE" reading of correlation plot viewer should be "correlation"
    And the "text of cell FLAT x AGE" reading of correlation plot viewer should be ""
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And the "error" reading of correlation plot viewer should be ""
    And no errors should have been logged
    When user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    And user removes "FLAT" column
    Then the "cells" reading of correlation plot viewer should be 16
    And no errors should have been logged

  @known-failure
  Scenario: Cells six times apart in coefficient are painted the same full red (GROK — color_coding.dart:322-337)
    # `_refreshGridColumns` sets minScale = -1 / maxScale = 1 on every matrix column
    # (correlation_plot_core.dart:195-196) so that a colour means the same thing everywhere, but
    # `ColorCoding.getGridCellAutoColor` takes the `isNumerical` branch (color_coding.dart:317-329)
    # and scales each grid column over its OWN min and max — the minScale/maxScale arguments live in
    # the branch below it, which a numerical column never reaches. So each column is stretched over
    # its own four values, the diagonal's stored 0 included: AGE x WEIGHT (+0.065) and
    # HEIGHT x WEIGHT (+0.412) are both #ff0000 because each is the largest value in its own column,
    # and a near-neutral +0.065 is painted as the extreme of a -1..1 scale. This is the assertion the
    # whole of correlation-plot-spec.ts was marked `test.fail(true, ...)` for, at line 343.
    Then the "correlation of WEIGHT and AGE" reading of correlation plot viewer should be between 0.0647 and 0.0649
    And the "correlation of WEIGHT and HEIGHT" reading of correlation plot viewer should be between 0.4124 and 0.4125
    And the "color of cell AGE x WEIGHT" and "color of cell HEIGHT x WEIGHT" readings of correlation plot viewer should differ
    And the "color of cell AGE x WEIGHT" reading of correlation plot viewer should not be "#ff0000"
    And no errors should have been logged
