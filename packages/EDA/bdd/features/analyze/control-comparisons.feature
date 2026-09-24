@eda @realizes:ml.menu.analyze.group-comparison.control-comparisons
Feature: Control comparisons
  ML | Analyze | Group Comparison | Control Comparisons... over demog: AGE of every RACE against the
  Asian control, Dunnett by default. Translated from the package's
  playwright/control-comparisons.test.ts; there is no Test Track case for it.

  The group sizes count nonmissing AGE values: 157 Black, 5266 Caucasian and 354 Other. What a
  filtered run should count is filtered-group-comparison.feature.

  Background:
    Given user is logged in
    And user opens demog dataset

  Scenario: Running the default comparison docks a box plot and the comparison statistics
    When user picks "ML > Analyze > Group Comparison > Control Comparisons..." from the top menu
    Then "Control comparisons" dialog should be visible
    And editor of Category input in "Control comparisons" dialog should have text "RACE"
    And Control input in "Control comparisons" dialog should have value "Asian"
    And editor of Feature input in "Control comparisons" dialog should have text "AGE"
    And Alpha input in "Control comparisons" dialog should have value "0.05"
    And Run button in "Control comparisons" dialog should be enabled
    When user clicks on Run button in "Control comparisons" dialog
    Then "Control comparisons" dialog should be hidden
    And box plot viewer should be visible
    And description of box plot viewer should be visible
    And description of box plot viewer should contain text "Asian"
    And box plot viewer should be painted
    And table "Control comparisons result" should be open
    And table "Control comparisons result" should have 3 rows
    And table "Control comparisons result" should have columns "Conclusion, Group, n, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g"
    And table "Control comparisons result" should have no missing values in "Mean diff" column
    And table "Control comparisons result" should have no missing values in "95% CI low" column
    And table "Control comparisons result" should have no missing values in "95% CI high" column
    And table "Control comparisons result" should have no missing values in "t" column
    And table "Control comparisons result" should have no missing values in "df" column
    And table "Control comparisons result" should have no missing values in "p (raw)" column
    And table "Control comparisons result" should have no missing values in "p (adj)" column
    And table "Control comparisons result" should have no missing values in "Hedges' g" column
    And second grid viewer should be bound to table "Control comparisons result"
    And the "text of cell 1 of Group" reading of second grid viewer should be "Black"
    And the "text of cell 1 of n" reading of second grid viewer should be "157"
    And the "text of cell 2 of Group" reading of second grid viewer should be "Caucasian"
    And the "text of cell 2 of n" reading of second grid viewer should be "5266"
    And the "text of cell 3 of Group" reading of second grid viewer should be "Other"
    And the "text of cell 3 of n" reading of second grid viewer should be "354"
    And no error or warning balloon should have been shown
    And no errors should have been logged
