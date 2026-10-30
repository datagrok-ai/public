@eda
Feature: Two-sample t-test compares a feature between the two groups of a category
  ML | Analyze | Group Comparison | T-test... opens the Two-sample t-test dialog on the current table
  (help `explore/group-comparison.md#t-test`): Category is the first column with exactly two groups,
  Feature the first number column, Alpha 0.05, Method Welch, Full report on. Run docks a box plot of the
  feature by the category with the conclusion and the p-value in its description, and, with Full report
  on, a one-row "T-test result" table under it. A category with three groups disables Run.

  Every scenario runs on a six-row table written here: group a holds v = 1, 2, 3 and group b 11, 12, 13
  (clearly different), w overlaps between the groups, site has three groups. Not a journey: every Run
  docks a box plot and a grid that the next scenario's viewer claims would resolve to. Left out on purpose:
  - the numbers of the test (t, df, p, CI, effect sizes) for Welch and Student: EDA src/tests
    (ttest-tests) are their place; Student is run once, without claiming its numbers;
  - ANOVA and Control Comparisons, the other items of the menu, have their own features
    (anova.feature, control-comparisons.feature);
  - the reason Run gives when it is disabled: it lives only in a tooltip (MISSING.md, needs a step).
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "groups" with:
      | grp | site | v  | w |
      | a   | x    | 1  | 5 |
      | a   | y    | 2  | 7 |
      | a   | z    | 3  | 6 |
      | b   | x    | 11 | 6 |
      | b   | y    | 12 | 5 |
      | b   | z    | 13 | 7 |

  Scenario: The dialog opens on the two-group column, and Run docks the box plot and the result table
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    Then "Two-sample t-test" dialog should be visible
    And editor of Category input in "Two-sample t-test" dialog should have text "grp"
    And editor of Feature input in "Two-sample t-test" dialog should have text "v"
    And Alpha input in "Two-sample t-test" dialog should have value "0.05"
    And "Full report" input in "Two-sample t-test" dialog should be checked
    And Run button in "Two-sample t-test" dialog should be enabled
    When user clicks on Run button in "Two-sample t-test" dialog
    Then "Two-sample t-test" dialog should be hidden
    And box plot viewer should be visible
    And description of box plot viewer should contain text "\"v\" differs between \"a\" and \"b\""
    And box plot viewer should be painted
    And table "T-test result" should be open
    And table "T-test result" should have 1 row
    And table "T-test result" should have columns "Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g"
    And second grid viewer should be bound to table "T-test result"
    And the "text of cell 1 of Conclusion" reading of second grid viewer should be "Significant"
    And table "groups" should have 6 rows
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Overlapping groups are reported as not differing
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    And user selects "w" in Feature input in "Two-sample t-test" dialog
    And user clicks on Run button in "Two-sample t-test" dialog
    Then box plot viewer should be visible
    And description of box plot viewer should contain text "\"w\" doesn't differ between \"a\" and \"b\""
    And the "text of cell 1 of Conclusion" reading of second grid viewer should be "Not significant"
    And no errors should have been logged

  Scenario: A category with three groups disables Run, and a two-group one enables it again
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    And user selects "site" in Category input in "Two-sample t-test" dialog
    Then Run button in "Two-sample t-test" dialog should be disabled
    When user selects "grp" in Category input in "Two-sample t-test" dialog
    Then Run button in "Two-sample t-test" dialog should be enabled
    When user clicks on CANCEL button in "Two-sample t-test" dialog
    Then "Two-sample t-test" dialog should be hidden
    And no errors should have been logged

  Scenario: Student with Full report off docks only the box plot
    When user picks "ML > Analyze > Group Comparison > T-test..." from the top menu
    And user selects "Student" in Method input in "Two-sample t-test" dialog
    And user unchecks "Full report" input in "Two-sample t-test" dialog
    And user clicks on Run button in "Two-sample t-test" dialog
    Then box plot viewer should be visible
    And description of box plot viewer should contain text "\"v\" differs between \"a\" and \"b\""
    And the open tableview should have 1 grid viewer
    And no errors should have been logged
