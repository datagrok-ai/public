@viewers @edge-cases @realizes:viewers.scatter-plot
Feature: NX edge cases — a column a calculated column reads is removed while a viewer shows the result
  The NX chain renames the columns its calculated columns read (`nx-chain.feature`); removing one
  is a rule of its own. When a column a formula reads is removed, the formula is recalculated in
  full and comes out empty, and the client warns with a balloon naming the dependent column
  (`core/docs/CALCULATED_COLUMNS.md`, "columns removed"); the column is removed through Remove in
  its header's context menu, the command that raises the balloon. A scatter plot that shows the
  calculated column stays open without an error. demog-1000; `AGE2` is `${AGE} * 2`, 52 in row 1 (AGE 26).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: Removing the column a calculated column reads empties the result and names it in a balloon
    When user adds a calculated column "AGE2" with formula "${AGE} * 2"
    Then the value of "AGE2" column in row 1 should be "52"
    And "AGE2" column should have no missing values
    When user adds a scatter plot viewer with:
      | xColumnName | AGE2   |
      | yColumnName | WEIGHT |
    And user picks "Remove" from the context menu of the "header AGE" area of grid
    Then the table should not have a column "AGE"
    And a warning balloon containing "AGE2 uses the removed column and lost its values" should have been shown
    And the value of "AGE2" column in row 1 should be ""
    And "AGE2" column should have missing values
    And scatter plot viewer should report no error
    And no errors should have been logged
