@journey @eda @realizes:ml.menu.analyze.group-comparison.control-comparisons
Feature: Group comparison of a filtered table
  A group comparison run while a filter is on should compare the rows that pass the filter.
  Translated from the GROK-20795 test of the package's playwright/control-comparisons.test.ts,
  which asserted the opposite on purpose so that the fix would break it.

  demog-1000 has 553 women (7 Asian, 18 Black, 480 Caucasian, 48 Other). Fixed on 2026-09-21: the dialog clones the category and feature columns
  through the table's filter before factorization. The group sizes it compares (the women only) are
  the EDA test "Control comparisons: dialog | A filtered table compares only the rows that pass the filter".

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: Running control comparisons preserves the filter and produces the comparison table
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    When user picks "ML > Analyze > Group Comparison > Control Comparisons..." from the top menu
    Then "Control comparisons" dialog should be visible
    And 553 rows should pass the filter
    When user clicks on Run button in "Control comparisons" dialog
    Then table "Control comparisons result" should be open
    And "Control comparisons" dialog should be hidden
    And 553 rows should pass the filter
    And table "Control comparisons result" should have 3 rows
    And table "Control comparisons result" should have no missing values in "n" column
    And second grid viewer should be bound to table "Control comparisons result"
    And the "text of cell 1 of Group" reading of second grid viewer should be "Black"
    And the "text of cell 2 of Group" reading of second grid viewer should be "Caucasian"
    And the "text of cell 3 of Group" reading of second grid viewer should be "Other"
    And no error or warning balloon should have been shown
    And no errors should have been logged
