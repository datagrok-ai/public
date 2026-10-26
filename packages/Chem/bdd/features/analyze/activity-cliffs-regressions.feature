@journey @realizes:chem.cp.activity-cliffs
Feature: Activity Cliffs regressions from closed issues
  The Activity cliffs panel the cliff link opens belongs to the view of its analysis: on activity-cliffs
  (29 molecules with an Activity column) Chem | Analyze | Activity Cliffs... with its defaults plots the
  cliffs, the link opens the panel, another table's view does not show it, and back on the analysed
  table's view it is there again (public#1512, reported on Bio's Sequence Activity Cliffs, which shares
  the panel).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens activity-cliffs dataset

  Scenario: The Activity cliffs panel stays with the view of its analysis (public#1512)
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    And user clicks on OK button in "Activity Cliffs" dialog
    Then a new column matching "^sali_" should have been added
    And the "cliffs" reading of scatter plot viewer should be 2
    When user clicks on "cliffs" button
    Then "Activity cliffs" dock panel should be visible
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    And "Activity cliffs" dock panel should not be visible
    Given user switches to the "activity_cliffs_test" table view
    Then the "activity_cliffs_test" view should be current
    And "Activity cliffs" dock panel should be visible
    And no errors should have been logged
