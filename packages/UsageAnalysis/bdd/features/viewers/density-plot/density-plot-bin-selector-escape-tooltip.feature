@journey @viewers @realizes:viewers.density-plot
Feature: Density plot bin selector, Escape and the bin tooltip
  Show Bin Selector puts the bin-count selector into the plot's layout and takes it out, as the
  layout resolved it (the reading is the layout's flag, not a check that the selector is on screen —
  see MISSING.md); Escape clears a selection made by clicking a bin; and the densest bin's tooltip
  gives its range of both binned columns (AGE 47..47, WEIGHT 62.20..64.40) and the 5 rows it holds.
  A column name with markup characters (`<b>X</b>`, AGE renamed) is shown in the bin tooltip as
  typed: the caption was built as HTML from the raw column names before core#1048.
  One journey on demog-1000 with a density plot of AGE against WEIGHT, as in
  density-plot-select-and-zoom.feature; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then the "rows shown" reading of density plot viewer should be 1000
    And density plot viewer should have a "densest bin" area

  Scenario: Show Bin Selector takes the bin selector out of the layout and puts it back, as the layout resolved it
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

  Scenario: A bin's tooltip gives the bin's range of both binned columns and its row count
    Then the "rows in densest bin" reading of density plot viewer should be 5
    When user hovers over the "densest bin" area of density plot viewer
    Then tooltip should be visible
    And tooltip should contain text "AGE: 47..47"
    And tooltip should contain text "WEIGHT: 62.20..64.40"
    And tooltip should contain text "5 rows"
    When user moves the pointer away from density plot viewer
    Then no errors should have been logged

  @realizes:core#1048
  Scenario: A column name with markup characters is shown in the bin tooltip as typed (core#1048)
    When user renames "AGE" column to "<b>X</b>"
    And user sets "xColumnName" property of density plot viewer to "<b>X</b>"
    And user hovers over the "densest bin" area of density plot viewer
    Then tooltip should contain text "WEIGHT"
    And tooltip should contain text "<b>X</b>"
    When user moves the pointer away from density plot viewer
    Then no errors should have been logged
    When user renames "<b>X</b>" column to "AGE"
    And user sets "xColumnName" property of density plot viewer to "AGE"
    Then the table should have a column "AGE"
    And the "rows in densest bin" reading of density plot viewer should be 5
    And no errors should have been logged
