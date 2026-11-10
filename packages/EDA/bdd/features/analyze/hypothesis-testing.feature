@eda
Feature: Two-sample t-test compares a feature between the two groups of a category
  ML | Analyze | Group Comparison | T-test... opens the Two-sample t-test dialog on the current table (help
  `explore/group-comparison.md#t-test`): Category is the first column with exactly two groups, Feature the first
  number column. A category with three groups disables Run. Run docks a box plot of the feature by the category with
  the conclusion in its description and, with Full report on (the default), a one-row "T-test result" table under it.

  Runs on a six-row table written here: group a holds v = 1, 2, 3 and group b 11, 12, 13; site has three groups.
  The test's numbers (Welch and Student) are EDA's ttest-tests; Full report off and a run on overlapping groups are
  one rule each (MISSING.md). Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "groups" with:
      | grp | site | v  |
      | a   | x    | 1  |
      | a   | y    | 2  |
      | a   | z    | 3  |
      | b   | x    | 11 |
      | b   | y    | 12 |
      | b   | z    | 13 |

  Scenario: The dialog opens on the two-group column, a three-group one disables Run, and Run docks the box plot and the result table
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    Then "Two-sample t-test" dialog should be visible
    And editor of Category input in "Two-sample t-test" dialog should have text "grp"
    And editor of Feature input in "Two-sample t-test" dialog should have text "v"
    When user selects "site" in Category input in "Two-sample t-test" dialog
    Then Run button in "Two-sample t-test" dialog should be disabled
    When user selects "grp" in Category input in "Two-sample t-test" dialog
    Then Run button in "Two-sample t-test" dialog should be enabled
    When user clicks on Run button in "Two-sample t-test" dialog
    Then description of box plot viewer should contain text "\"v\" differs between \"a\" and \"b\""
    And table "T-test result" should have 1 row
    And table "T-test result" should have columns "Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g"
    And second grid viewer should be bound to table "T-test result"
    And no error or warning balloon should have been shown
    And no errors should have been logged
