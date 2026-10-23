@journey @viewers @realizes:viewers.density-plot
Feature: Density plot bin selector, Escape and the bin tooltip
  Show Bin Selector puts the bin-count selector on the plot and takes it away, as the layout
  resolved it; Escape clears a selection made by clicking a bin; and a bin's tooltip names the range
  of both binned columns. One journey on demog-1000 with a density plot of AGE against WEIGHT, as in
  density-plot-select-and-zoom.feature; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then the "rows shown" reading of density plot viewer should be 1000
    And density plot viewer should have a "densest bin" area

  Scenario: Show Bin Selector takes the bin selector off the plot and puts it back
    Then the "bin selector shown" reading of density plot viewer should be "true"
    When user sets "showBinSelector" property of density plot viewer to "false"
    Then the "bin selector shown" reading of density plot viewer should be "false"
    And the "x selector shown" reading of density plot viewer should be "true"
    When user sets "showBinSelector" property of density plot viewer to "true"
    Then the "bin selector shown" reading of density plot viewer should be "true"
    And no errors should have been logged

  Scenario: Escape clears the selection a bin click made
    Given user clears the row selection
    When user clicks on the "densest bin" area of density plot viewer
    Then the "rows selected" and "rows in densest bin" readings of density plot viewer should be the same
    And some rows should be selected
    When user presses Escape
    Then no rows should be selected
    And the "rows selected" reading of density plot viewer should be 0
    And no errors should have been logged

  Scenario: A bin's tooltip names the range of both binned columns
    When user hovers over the "densest bin" area of density plot viewer
    Then tooltip should be visible
    And tooltip should contain text "AGE"
    And tooltip should contain text "WEIGHT"
    And tooltip should contain text "rows"
    When user moves the pointer away from grid
    Then no errors should have been logged
