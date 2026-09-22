@journey @realizes:chem.cp.r-group-analysis
Feature: R-Groups Analysis with MCS, Replace latest and no core
  On five copies of benzene, R-Groups Analysis with the MCS core and Visual analysis on adds nothing,
  and on the 1000 unrelated molecules of smiles, once MCS has drawn its core, it says "No R-Groups
  were found". On sar-small the same run adds R-group columns and a trellis plot; run
  again with Replace latest off it adds a second set beside the first, and run again with Replace
  latest on it takes the latest set away and puts its own in place. A run with an empty sketcher
  says "No core was provided" and keeps the columns of the run before it.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: MCS over molecules that are all the same finds no R-groups
    Given user opens a table "same_mols" with:
      | id | smiles   |
      | 1  | c1ccccc1 |
      | 2  | c1ccccc1 |
      | 3  | c1ccccc1 |
      | 4  | c1ccccc1 |
      | 5  | c1ccccc1 |
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "R-Groups Analysis" dialog should be visible
    And "Visual analysis" input in "R-Groups Analysis" dialog should be checked
    When user clicks on MCS button in "R-Groups Analysis" dialog
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then the top menu command should have completed
    And no new column should have been added
    And the table should have 5 rows
    And no errors should have been logged

  Scenario: MCS over one series adds R-group columns and a trellis plot
    Given user opens sar-small dataset
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then the top menu command should have completed
    And a new column matching "^R1" should have been added
    And trellis plot viewer should be visible
    And no errors should have been logged

  Scenario: With Replace latest off the second run adds its own columns beside the first
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "Replace latest" input in "R-Groups Analysis" dialog should be visible
    When user clicks on MCS button in "R-Groups Analysis" dialog
    And user unchecks "Replace latest" input in "R-Groups Analysis" dialog
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then the top menu command should have completed
    And a new column matching "^R1" should have been added
    And the open tableview should have 2 trellis plot viewers
    And no errors should have been logged

  Scenario: With Replace latest on the third run takes the latest set away
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    And user checks "Replace latest" input in "R-Groups Analysis" dialog
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then the top menu command should have completed
    And the open tableview should have 2 trellis plot viewers
    And no errors should have been logged


  Scenario: A run without a core says so and keeps the results of the run before it
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No core was provided" should have been shown
    And the table should have a column "R1"

  Scenario: MCS over a table of unrelated molecules finds no R-groups
    Given user opens smiles dataset
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    Then the canvases of "R-Groups Analysis" dialog should be painted in at least 2 colors
    When user clicks on OK button in "R-Groups Analysis" dialog
    Then the top menu command should have completed
    And an error balloon containing "No R-Groups were found" should have been shown
    And no new column should have been added
    And the table should have 1000 rows
