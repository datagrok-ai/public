@journey @realizes:chem.cp.mmp-analysis
Feature: Matched Molecular Pairs over molecules and one activity
  Chem | Analyze | Matched Molecular Pairs... on sar-small (200 molecules of one series) opens on the
  smiles column; LD(50) is picked as its activity. OK builds the Matched Molecular Pairs Analysis
  viewer: it names the activity it ran on, holds substitutions and molecule pairs, offers the four
  tabs, and fills the Generation tab.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens sar-small dataset

  Scenario: The dialog opens on the molecules column and the activity
    When user picks "Chem > Analyze > Matched Molecular Pairs..." from the top menu
    Then "Matched Molecular Pairs" dialog should be visible
    And Column input in "Matched Molecular Pairs" dialog should contain text "smiles"
    And Activities input in "Matched Molecular Pairs" dialog should contain text "Activities(0)"
    When user clicks on editor of Activities input in "Matched Molecular Pairs" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    Then the "text of cell 6 of __name" reading of grid viewer in "Select columns..." dialog should be "LD(50)"
    When user clicks on OK button in "Select columns..." dialog
    Then Activities input in "Matched Molecular Pairs" dialog should contain text "Activities(1)"
    And Scaling input in "Matched Molecular Pairs" dialog should be visible
    And no errors should have been logged

  Scenario: The run builds the viewer with its substitutions and pairs
    When user clicks on OK button in "Matched Molecular Pairs" dialog
    Then Matched Molecular Pairs Analysis viewer should be visible
    And Matched Molecular Pairs Analysis viewer should have finished its analysis
    And the "activities" reading of Matched Molecular Pairs Analysis viewer should be "LD(50)"
    And the "molecules column" reading of Matched Molecular Pairs Analysis viewer should be "smiles"
    And the "substitutions" reading of Matched Molecular Pairs Analysis viewer should be at least 1
    And the "pairs" reading of Matched Molecular Pairs Analysis viewer should be at least 1
    And the table should have 200 rows
    And no errors should have been logged

  Scenario: The tabs of the viewer are its four analyses
    Then Substitutions tab in Matched Molecular Pairs Analysis viewer should be visible
    And Fragments tab in Matched Molecular Pairs Analysis viewer should be visible
    And Cliffs tab in Matched Molecular Pairs Analysis viewer should be visible
    And Generation tab in Matched Molecular Pairs Analysis viewer should be visible
    And the "tab" reading of Matched Molecular Pairs Analysis viewer should be "Substitutions"
    And no errors should have been logged

  Scenario: The Generation tab fills its grid
    When user clicks on Generation tab in Matched Molecular Pairs Analysis viewer
    Then the "tab" reading of Matched Molecular Pairs Analysis viewer should be "Generation"
    And the "generated" reading of Matched Molecular Pairs Analysis viewer should be at least 1
    And no errors should have been logged
