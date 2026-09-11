@eda @realizes:ml.menu.analyze.group-comparison.control-comparisons
Feature: Control comparisons
  ML | Analyze | Group Comparison | Control Comparisons... over demog: AGE of every RACE against the
  Asian control, Dunnett by default. Translated from the package's
  playwright/control-comparisons.test.ts; there is no Test Track case for it.

  The group sizes are demog's own: 157 Black, 5266 Caucasian and 354 Other. What a filtered run
  should count is filtered-group-comparison.feature.

  Background:
    Given user is logged in
    And user opens demog dataset

  Scenario: The dialog opens on a category, its control and a feature
    When user picks "ML > Analyze > Group Comparison > Control Comparisons..." from the top menu
    Then "Control comparisons" dialog should be visible
    And editor of Category input in "Control comparisons" dialog should have text "RACE"
    And Control input in "Control comparisons" dialog should have value "Asian"
    And editor of Feature input in "Control comparisons" dialog should have text "AGE"
    And Alpha input in "Control comparisons" dialog should have value "0.05"

  Scenario: Running it docks a box plot and the table of the comparisons
    When user picks "ML > Analyze > Group Comparison > Control Comparisons..." from the top menu
    And user clicks on Run button in "Control comparisons" dialog
    Then the top menu command should have completed
    And "Control comparisons" dialog should be hidden
    And box plot viewer should be visible
    And "Description" property of box plot viewer should contain "Asian"
    And box plot viewer should be painted
    And table "Control comparisons result" should be open
    And table "Control comparisons result" should have 3 rows
    And table "Control comparisons result" should have columns "Conclusion, Group, n, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g"
    And second grid viewer should be bound to table "Control comparisons result"
    And the "text of cell 1 of Group" reading of second grid viewer should be "Black"
    And the "text of cell 1 of n" reading of second grid viewer should be "157"
    And the "text of cell 2 of n" reading of second grid viewer should be "5266"
    And the "text of cell 3 of n" reading of second grid viewer should be "354"
    And no error or warning balloon should have been shown
    And no errors should have been logged
