@journey @realizes:chem.int.empty-input-analyses
Feature: R-Groups Analysis on an all-empty molecule column
  A ten-row table whose structure column is a Molecule column with nothing in it. R-Groups Analysis
  with the MCS strategy derives no core from it, says so in a balloon and adds no column, and the
  grid keeps working afterwards (GROK-16329).

  Chemical Space on such a column is tested in Chem src/tests/menu-tests-chem-space.ts 'chemSpace.allEmpty'.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens a table "empty_mols" with:
      | id | structure |
      | 1  |           |
      | 2  |           |
      | 3  |           |
      | 4  |           |
      | 5  |           |
      | 6  |           |
      | 7  |           |
      | 8  |           |
      | 9  |           |
      | 10 |           |
    And user sets the semantic type of "structure" column to "Molecule"

  Scenario: R-Groups Analysis with MCS says it has no core and adds nothing
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "R-Groups Analysis" dialog should be visible
    When user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No core was provided" should have been shown
    And no new column should have been added
    And the table should have 10 rows
    And all rows should pass the filter
    When user clicks on the "cell 4 of structure" area of grid
    Then row 4 should be current
    And the current column should be "structure"
    And no errors should have been logged
