@viewers @realizes:viewers.filters
Feature: Expression filter over All Columns
  With Column set to All Columns, the expression card looks for its value in every column it can
  compare it with, and keeps the rows where any of them holds it. The table is written below, and
  the Background checks the type each column gets from it: Name and Prediction are strings, Score
  and Mass doubles, Idea ID an integer, Registered a date. Score's short decimals are held in 32
  bits, so its first row reads 0.4699999988079071, while its 0.5 is exact. Prediction is Medium in
  2 rows (Alpha, Delta); Idea ID is 634783 in 2 (Alpha, Epsilon) and above 1000 in 3 (Alpha, Beta,
  Epsilon); Mass is 15.05460453 in 1 (Alpha); Score is 0.5 in 1 (Beta) and shows 0.47 in 2 (Alpha,
  Gamma); Registered is 2020-01-15 in 2 (Alpha, Delta). No value is held by two columns, so each
  search has one right answer.
  Every scenario claims the rule the card committed as well as the rows it keeps, since the card
  filters by the value while it is being typed. Text matches with "equals" and with "=", in any
  case, and "equals" ignores spaces around the value. A number matches with "=" whether it is an
  integer or a double, and a comparison reaches every numeric column and no date. A value no column
  holds keeps no row, and Remove All brings every row back. The date is claimed by its count and by
  its rows through Prediction, whose Medium rows are the same two; there is no row check for dates.
  "equals" with a number reaches the numeric columns in All Columns mode, and 0.47 with "=" keeps
  the two rows that show it, in All Columns mode and on Score itself: the card compares a typed
  number in the 32-bit form Score stores (GROK-20977).

  Background:
    Given user is logged in
    And user opens a table "expression-types" with:
      | Name    | Prediction | Score | Mass        | Idea ID | Registered |
      | Alpha   | Medium     | 0.47  | 15.05460453 | 634783  | 2020-01-15 |
      | Beta    | Low        | 0.5   | 2.274373293 | 634784  | 2021-06-30 |
      | Gamma   | High       | 0.47  | 10.34033203 | 12      | 2022-03-01 |
      | Delta   | Medium     | 0.83  | 12.90704918 | 7       | 2020-01-15 |
      | Epsilon | Low        | 1.25  | 1.877748251 | 634783  | 2023-11-20 |
    Then "Prediction" column should have type "string"
    And "Score" column should have type "double"
    And "Mass" column should have type "double"
    And "Idea ID" column should have type "int"
    And "Registered" column should have type "datetime"
    And the value of "Score" column in row 1 should be "0.4699999988079071"
    And the value of "Score" column in row 2 should be "0.5"
    When user opens an empty filter panel
    And user picks "Add Filter | Expression" from the filter panel menu
    Then "Expression" filter card should be visible
    And all rows should pass the filter

  Scenario Outline: A text value keeps the rows of its category, whatever the operation or case
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "<operation>" in Operation input in "Expression" filter card
    And user types "<value>" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  <operation> <value>"
    And 2 rows should pass the filter
    And the filter should pass exactly the rows where "Prediction" is "Medium"
    And no errors should have been logged

    Examples:
      | operation | value  |
      | equals    | Medium |
      | =         | Medium |
      | equals    | medium |
      | =         | MEDIUM |

  Scenario: Spaces around a text value are not part of it
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types " Medium " into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  equals  Medium "
    And 2 rows should pass the filter
    And the filter should pass exactly the rows where "Prediction" is "Medium"
    And no errors should have been logged

  Scenario Outline: A number with "=" keeps the rows that hold it, integer or double
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "=" in Operation input in "Expression" filter card
    And user types "<value>" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  = <value>"
    And <count> rows should pass the filter
    And the filter should pass exactly the rows where "<column>" is "<value>"
    And no errors should have been logged

    Examples:
      | value       | column  | count |
      | 634783      | Idea ID | 2     |
      | 15.05460453 | Mass    | 1     |
      | 0.5         | Score   | 1     |

  Scenario: A comparison reaches every numeric column and no date
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects ">" in Operation input in "Expression" filter card
    And user types "1000" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  > 1000"
    And 3 rows should pass the filter
    And the filter should pass exactly the rows where "Idea ID" is between 1001 and 1000000
    And no errors should have been logged

  Scenario: A date keeps the rows registered on it
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "2020-01-15" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  equals 2020-01-15"
    And 2 rows should pass the filter
    And the filter should pass exactly the rows where "Prediction" is "Medium"
    And no errors should have been logged

  Scenario: A value no column holds keeps no row, and Remove All brings every row back
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "Nothing" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  equals Nothing"
    And 0 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario Outline: "equals" with a number keeps the rows that hold it
    When user selects "All Columns" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "<value>" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "*  equals <value>"
    And <count> rows should pass the filter
    And the filter should pass exactly the rows where "<column>" is "<value>"

    Examples:
      | value       | column  | count |
      | 634783      | Idea ID | 2     |
      | 15.05460453 | Mass    | 1     |

  Scenario Outline: A short decimal with "=" keeps the rows that show it
    When user selects "<column mode>" in Column input in "Expression" filter card
    And user selects "=" in Operation input in "Expression" filter card
    And user types "0.47" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "<rule>"
    And 2 rows should pass the filter
    And the filter should pass exactly the rows where "Score" is between 0.465 and 0.475

    Examples:
      | column mode | rule            |
      | All Columns | *  = 0.47       |
      | Score       | ${Score} = 0.47 |
