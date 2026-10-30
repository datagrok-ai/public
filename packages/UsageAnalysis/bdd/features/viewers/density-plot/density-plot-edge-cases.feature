@viewers @realizes:viewers.density-plot
Feature: Density plot edge cases
  Edge cases of the density plot that the density plot audit core#1048 (open) found wrong, each
  written as the behaviour the audit and the plot's own contract expect and tagged as a known failure
  of core#1048: a click on a rectangular bin selects the rows counted in it (as it does for hexagons,
  `density-plot-select-and-zoom.feature`); a change of Bins after a zoom bins the same rows while
  Bin To Range is off (the column's range is binned unless Bin To Range is on); a logarithmic X axis
  leaves out the rows it cannot place, so only the positive values are binned; and a column name with
  markup characters is shown in the bin tooltip as typed.

  Background:
    Given user is logged in

  @edge-cases @realizes:core#1048
  Scenario: A click on a rectangular bin selects exactly the rows counted in it
    Guards core#1048 ("Wrong 1", fixed): with rectangular bins the click selected a different number
    of rows than the bin held.
    Given user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE       |
      | yColumnName | WEIGHT    |
      | binShape    | rectangle |
    Then the "bin shape" reading of density plot viewer should be "rectangle"
    And density plot viewer should have a "densest bin" area
    When user clears the row selection
    And user clicks on the "densest bin" area of density plot viewer
    Then the "rows in densest bin" reading of density plot viewer should be at least 1
    And the "rows selected" and "rows in densest bin" readings of density plot viewer should be the same

  @edge-cases @realizes:core#1048
  Scenario: Changing Bins after a zoom bins the same rows while Bin To Range is off
    Guards core#1048 ("Wrong 2", fixed): with rectangular bins a change of Bins after a zoom dropped
    rows from the plot.
    Given user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE       |
      | yColumnName | WEIGHT    |
      | binShape    | rectangle |
    Then the "rows shown" reading of density plot viewer should be 1000
    And the "bin to range" reading of density plot viewer should be "false"
    When user drags a zoom box over the "view" area of density plot viewer
    And user sets "bins" property of density plot viewer to "30"
    Then the "rows shown" reading of density plot viewer should be 1000

  @edge-cases @realizes:core#1048
  Scenario: A logarithmic X axis bins only the positive values
    Guards core#1048 ("Wrong 5", fixed): with hexagons the rows with zero or negative X were binned at
    the bottom of the log axis instead of being left out.

    Given user opens a table "log-x" with:
      | X   | Y |
      | -5  | 1 |
      | -1  | 2 |
      | 0   | 3 |
      | 2   | 4 |
      | 10  | 5 |
      | 100 | 6 |
    And user adds a density plot viewer with:
      | xColumnName | X |
      | yColumnName | Y |
    Then the "rows shown" reading of density plot viewer should be 6
    When user sets "xAxisType" property of density plot viewer to "logarithmic"
    Then the "rows shown" reading of density plot viewer should be 3

  @edge-cases @realizes:core#1048
  Scenario: A column name with markup characters is shown in the bin tooltip as typed
    Guards core#1048 ("Other", fixed): the tooltip caption was built as HTML from the raw column
    names, so "<b>X</b>" showed as a bold X.
    Given user opens demog-1000 dataset
    And user renames "AGE" column to "<b>X</b>"
    And user adds a density plot viewer with:
      | xColumnName | <b>X</b> |
      | yColumnName | WEIGHT   |
    When user hovers over the "densest bin" area of density plot viewer
    Then tooltip should contain text "WEIGHT"
    And tooltip should contain text "<b>X</b>"
