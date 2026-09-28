@journey @realizes:chem.cp.r-group-analysis
Feature: R-Groups Analysis with MCS, Replace latest and no core
  On five copies of benzene, R-Groups Analysis with the MCS core and Visual analysis on says "No
  R-Groups were found" and adds no R-group column (the decomposition's Core, highlight and isMatch
  columns come all the same), and on the 1000 unrelated molecules of smiles, once MCS has drawn its
  core, it says the same and adds nothing. On sar-small the same run adds R-group columns and a trellis plot;
  run again with Replace latest off it adds a second set (R1_1, ...) beside the first, and run again
  with Replace latest on it takes the latest set away and puts its own in place under the same names.
  A run with an empty sketcher says "No core was provided" and keeps the latest set and its plot.

  The menu command only opens the dialog, so no claim waits on it: each run is over when its result
  is there — the balloon, the columns, or the task bar's "R-Group analysis running..." entry gone.

  The dialog remembers Only match at R groups per account, and an MCS core carries no labelled R
  group, so the first run unticks it through the dialog's gear; the account is left with the option
  off, the product's default.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
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
    # the dialog remembers Only match at R groups per account, and an MCS core carries no labelled R group
    When user clicks on R-Groups settings icon
    And user unchecks "Only match at R groups" input in "R-Groups Analysis" dialog
    When user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No R-Groups were found" should have been shown
    And the table should not have a column "R1"
    And the table should have 5 rows
    And no errors should have been logged

  Scenario: MCS over one series adds R-group columns and a trellis plot
    Given user opens sar-small dataset
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then a new column matching "^R1" should have been added
    And trellis plot viewer should be visible
    And no errors should have been logged

  Scenario: With Replace latest off the second run adds its own columns beside the first
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "Replace latest" input in "R-Groups Analysis" dialog should be visible
    When user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    And user unchecks "Replace latest" input in "R-Groups Analysis" dialog
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then a new column "R1_1" should have been added
    And the open tableview should have 2 trellis plot viewers
    And no errors should have been logged

  Scenario: With Replace latest on the third run takes the latest set away
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    And user checks "Replace latest" input in "R-Groups Analysis" dialog
    And user watches the task bar
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then the task bar should have finished "R-Group analysis running"
    And the open tableview should have 2 trellis plot viewers
    And the table should have a column "R1_1"
    And no new column should have been added
    And no errors should have been logged

  Scenario: A run without a core says so and keeps the latest results
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No core was provided" should have been shown
    And the table should have a column "R1_1"
    And the open tableview should have 2 trellis plot viewers

  Scenario: MCS over a table of unrelated molecules finds no R-groups
    Given user opens smiles dataset
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    And user clicks on MCS button in "R-Groups Analysis" dialog
    And "R-Groups Analysis" dialog should have finished updating
    Then the canvases of "R-Groups Analysis" dialog should be painted in at least 2 colors
    When user clicks on OK button in "R-Groups Analysis" dialog
    Then an error balloon containing "No R-Groups were found" should have been shown
    And no new column should have been added
    And the table should have 1000 rows
