@journey @eda @realizes:ml.menu.analyze.group-comparison.control-comparisons
Feature: Group comparison of a filtered table
  A group comparison run while a filter is on should compare the rows that pass the filter.
  Translated from the GROK-20795 test of the package's playwright/control-comparisons.test.ts,
  which asserted the opposite on purpose so that the fix would break it.

  demog has 3243 women: 104 Black, 2823 Caucasian and 279 Other, against 157, 5266 and 354 on the
  whole table. The dialog counts the whole table today, so the scenario is @known-failure; when it
  passes, GROK-20795 is fixed and the tag goes.

  Background:
    Given user is logged in
    And user opens demog dataset

  @known-failure
  Scenario: Control comparisons count only the women when the table is filtered to them
    When user filters rows where "SEX" is "F"
    Then 3243 rows should pass the filter
    When user picks "ML > Analyze > Group Comparison > Control Comparisons..." from the top menu
    Then 3243 rows should pass the filter
    When user clicks on Run button in "Control comparisons" dialog
    Then the top menu command should have completed
    And 3243 rows should pass the filter
    And second grid viewer should be bound to table "Control comparisons result"
    And the "text of cell 1 of Group" reading of second grid viewer should be "Black"
    And the "text of cell 1 of n" reading of second grid viewer should be "104"
    And the "text of cell 2 of n" reading of second grid viewer should be "2823"
    And the "text of cell 3 of n" reading of second grid viewer should be "279"
