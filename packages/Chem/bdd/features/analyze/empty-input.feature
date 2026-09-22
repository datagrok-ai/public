@journey @realizes:chem.int.empty-input-analyses
Feature: R-Groups Analysis and Chemical Space on an all-empty molecule column
  A ten-row table whose structure column is a Molecule column with nothing in it. R-Groups Analysis
  with the MCS strategy derives no core from it, says so in a balloon and adds no column, and the
  grid keeps working afterwards (GROK-16329). Chemical Space runs on the same column: it adds the
  pair of embedding columns, a cluster column and the scatter plot of the embedding.

  Background:
    Given user is logged in
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
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No core was provided" should have been shown
    And no new column should have been added
    And the table should have 10 rows
    And all rows should pass the filter
    When user clicks on the "cell 4 of structure" area of grid
    Then row 4 should be current
    And the current column should be "structure"
    And no errors should have been logged

  Scenario: Chemical Space embeds the empty column and plots it
    When user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    And Column input in "Chem Space" dialog should contain text "structure"
    When user clicks on OK button in "Chem Space" dialog
    Then the "Chem Space" dialog should close
    And the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And a new column matching "^Embed_Y_" should have been added
    And the newest column matching "^Embed_X_" should have no missing values
    And a new column matching "^Cluster " should have been added
    And scatter plot viewer should be visible
    And the table should have 10 rows
    And no errors should have been logged
