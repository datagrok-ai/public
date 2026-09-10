@journey @viewers @realizes:viewers.correlation-plot
Feature: Correlation plot — which rows it correlates, which table, and what survives a round-trip
  The coefficient is not a property of two columns alone: it is computed over the rows the plot's
  Row Source and its own filter leave, and it moves when the table's filter moves. `rows shown` says
  how many those are and `correlation of X and Y` says what they produced, so every scenario here is
  a pair of numbers rather than "the canvas differs by more than the idle noise floor", which is what
  the spec this replaces measured — with the row-source case tuned to a `> 2000` pixel delta on a
  5850-row table and the viewer-filter case to `> 100`.
  The cross-check against `DG.Stats` follows the row source too: it clones the table over the same
  bitset the plot uses, so "the matrix agrees with an independent Pearson" stays a real claim when
  only twenty rows are selected.
  The font scenario ends the journey because the product does not undo it: `row height` starts at
  20, the field's initialiser, and becomes `parseSize(defaultCellFont) x 1.4` the first time the
  font is written — so putting the font back gives 18.2, not 20. The old spec re-ran its whole
  geometry calibration after this step for exactly that reason.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then 1000 rows should pass the filter
    And the "rows shown" reading of correlation plot viewer should be 1000
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And correlation plot viewer should be painted

  Scenario: The viewer's own filter narrows what it correlates and leaves the table alone
    When user sets "filter" property of correlation plot viewer to "${AGE} > 40"
    Then the "rows shown" reading of correlation plot viewer should be 635
    And 1000 rows should pass the filter
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2040 and -0.2039
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.20"
    And correlation plot viewer should have repainted
    When user sets "filter" property of correlation plot viewer to ""
    Then the "rows shown" reading of correlation plot viewer should be 1000
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And no errors should have been logged

  Scenario: A filter on the table moves the coefficients with it
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of correlation plot viewer should be 447
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should differ from before
    And the correlation of "HEIGHT" and "AGE" of correlation plot viewer should match the Pearson coefficient of the table
    And correlation plot viewer should have repainted
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of correlation plot viewer should be 1000
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And no errors should have been logged

  Scenario: Row Source Selected correlates the selected rows and nothing else
    When user selects the first 20 rows
    Then 20 rows should be selected
    When user sets "rowSource" property of correlation plot viewer to "Selected"
    Then the "rows shown" reading of correlation plot viewer should be 20
    And the correlation of "HEIGHT" and "AGE" of correlation plot viewer should match the Pearson coefficient of the table
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2161 and -0.2160
    And correlation plot viewer should have repainted
    When user sets "rowSource" property of correlation plot viewer to "Filtered"
    And user clears the row selection
    Then the "rows shown" reading of correlation plot viewer should be 1000
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And no errors should have been logged

  Scenario: Bound to another table the matrix becomes that table's numerical columns
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of correlation plot viewer to "spgi-100"
    Then correlation plot viewer should be bound to table "spgi-100"
    And the "rows shown" reading of correlation plot viewer should be 100
    And the "cells" reading of correlation plot viewer should be 529
    And the "numerical columns" reading of correlation plot viewer should contain "TPSA"
    And the "numerical columns" reading of correlation plot viewer should not contain "HEIGHT"
    And the "error" reading of correlation plot viewer should be ""
    And correlation plot viewer should be painted
    When user sets "table" property of correlation plot viewer to "demog-1000"
    Then correlation plot viewer should be bound to table "demog-1000"
    And the "cells" reading of correlation plot viewer should be 16
    And the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And no errors should have been logged

  Scenario: A saved layout brings the configured matrix back
    When user sets properties of correlation plot viewer:
      | correlationType | Spearman            |
      | showPearsonR    | false               |
      | xColumnNames    | AGE, HEIGHT, WEIGHT |
      | yColumnNames    | AGE, HEIGHT         |
    Then the "cells" reading of correlation plot viewer should be 6
    And user saves the layout of the current table view
    When user sets properties of correlation plot viewer:
      | correlationType | Pearson                      |
      | showPearsonR    | true                         |
      | xColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of correlation plot viewer should be 16
    When user loads the saved layout
    Then correlation plot viewer should be visible
    And the "cells" reading of correlation plot viewer should be 6
    And the "correlation type" reading of correlation plot viewer should be "Spearman"
    And the "show pearson r" reading of correlation plot viewer should be "false"
    And the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT"
    And the "y columns" reading of correlation plot viewer should be "AGE, HEIGHT"
    And the correlation of "WEIGHT" and "HEIGHT" of correlation plot viewer should match the Spearman coefficient of the table
    When user sets properties of correlation plot viewer:
      | correlationType | Pearson                      |
      | showPearsonR    | true                         |
      | xColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of correlation plot viewer should be 16
    And no errors should have been logged

  Scenario: A project round-trip brings it back too
    When user sets properties of correlation plot viewer:
      | correlationType | Spearman            |
      | xColumnNames    | AGE, HEIGHT, WEIGHT |
      | yColumnNames    | AGE, HEIGHT         |
    Then the "cells" reading of correlation plot viewer should be 6
    When user saves the current view as project "bdd correlation matrix"
    And user closes all views
    And user opens the "bdd correlation matrix" project
    Then correlation plot viewer should be visible
    And the "cells" reading of correlation plot viewer should be 6
    And the "correlation type" reading of correlation plot viewer should be "Spearman"
    And the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT"
    And the "rows shown" reading of correlation plot viewer should be 1000
    And the correlation of "WEIGHT" and "HEIGHT" of correlation plot viewer should match the Spearman coefficient of the table
    When user sets properties of correlation plot viewer:
      | correlationType | Pearson                      |
      | xColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames    | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of correlation plot viewer should be 16
    And no errors should have been logged

  Scenario: A bigger cell font makes the rows taller and moves no coefficient
    Then the "row height" reading of correlation plot viewer should be 20
    When user sets "defaultCellFont" property of correlation plot viewer to 'normal normal 20px "Roboto"'
    Then the "row height" reading of correlation plot viewer should be between 27.9 and 28.1
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And correlation plot viewer should have repainted
    When user sets "defaultCellFont" property of correlation plot viewer to 'normal normal 13px "Roboto"'
    Then the "row height" reading of correlation plot viewer should be between 18.1 and 18.3
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And no errors should have been logged
