@journey @viewers @realizes:viewers.filters
Feature: Numeric matching in the search box and the expression filter
  A float column keeps its values in single precision, so the cell that shows as 0.47 holds
  0.4699999988079071; the search box and the expression filter also compare a typed number in the
  form the column stores it, so "0.47" finds both 0.47 rows, ">= 0.83" keeps the 0.83 row and
  "< 0.47" drops it. The expression filter's All Columns mode offers the word "equals" next to the numeric
  "=", and both reach the numeric columns: on spgi, "equals 634783" keeps the one row whose
  CAST Idea ID is 634783 (its Id, "CAST-634783", is a string and does not equal the number).
  Not translated: which rows the search keeps is claimed by count only — the row check compares
  the stored single-precision value against a typed bound, which is the very comparison under test.

  Background:
    Given user is logged in
    And the toolbox pane is shown
    And user opens a table "scores" with:
      | Score | Prediction |
      | 0.47  | Medium     |
      | 0.5   | Low        |
      | 0.47  | High       |
      | 0.83  | Medium     |
      | 1.25  | Low        |
    Then all rows should pass the filter

  Scenario: Ctrl+F, a float typed as the cell shows it, Enter — the search box filters the rows
    When user presses Control+f in grid overlay
    Then table search should be visible
    When user types "0.47" into table search
    And user presses Enter in table search
    Then 2 rows should pass the filter
    When user types ">= 0.83" into table search
    And user presses Enter in table search
    Then 2 rows should pass the filter
    When user types "< 0.47" into table search
    And user presses Enter in table search
    Then 0 rows should pass the filter
    When user types "0.47-0.5" into table search
    And user presses Enter in table search
    Then 3 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: The expression filter compares a float column against the value it shows
    When user opens an empty filter panel
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "Score" in Column input in "Expression" filter card
    And user selects "=" in Operation input in "Expression" filter card
    And user types "0.47" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 2 rows should pass the filter
    And the "categories of Expression" reading of filter panel should be "${Score} = 0.47"
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "Score" in Column input in "Expression" filter card
    And user selects ">=" in Operation input in "Expression" filter card
    And user types "0.83" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 2 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "0.47" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 2 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: All Columns with the word "equals" reaches the numeric columns of spgi
    When user opens spgi dataset
    And user opens an empty filter panel
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "634783" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 1 row should pass the filter
    And the filter should pass exactly the rows where "CAST Idea ID" is between 634783 and 634783
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged
