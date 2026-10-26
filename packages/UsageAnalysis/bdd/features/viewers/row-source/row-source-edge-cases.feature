@viewers @edge-cases @realizes:viewers.scatter-plot @realizes:viewers.histogram @realizes:viewers.pie-chart
Feature: Row Source edge cases — a Filter formula that is not a condition, a row set nothing is left in
  A viewer's Filter formula that does not evaluate to true or false (`${AGE} + 1`) is refused with
  the "Invalid formula" error balloon, and the viewer draws its row source as if no formula were set
  (`DataFrameViewer.createCombinedRowFilter`); a formula that is a condition applied after it works
  again. Under Row Source FilteredSelected a viewer draws the rows that are both filtered and
  selected (`table-view-1.md`, Row Source), so with every selected row filtered out it draws
  none, and logs no error. demog-1000: 1000 rows, 519 with AGE > 44, 553 F and 447 M; AGE and
  WEIGHT have no blanks.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario Outline: A <viewer> Filter formula that is not a condition is refused with an error balloon and ignored
    When user adds a <viewer> viewer with:
      | <column property> | AGE |
    Then the "rows shown" reading of <viewer> viewer should be 1000
    When user sets "filter" property of <viewer> viewer to "${AGE} + 1"
    Then an error balloon containing "Invalid formula" should have been shown
    And the "rows shown" reading of <viewer> viewer should be 1000
    When user sets "filter" property of <viewer> viewer to "${AGE} > 44"
    Then the "rows shown" reading of <viewer> viewer should be 519
    And no errors should have been logged

    Examples:
      | viewer       | column property    |
      | scatter plot | xColumnName        |
      | histogram    | valueColumnName    |
      | pie chart    | categoryColumnName |

  Scenario: Under FilteredSelected with every selected row filtered out the viewer draws no row and logs no error
    When user adds a scatter plot viewer with:
      | xColumnName | AGE              |
      | yColumnName | WEIGHT           |
      | rowSource   | FilteredSelected |
    And user selects rows where "SEX" is "M"
    Then 447 rows should be selected
    And the "rows shown" reading of scatter plot viewer should be 447
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And the "rows shown" reading of scatter plot viewer should be 0
    And no errors should have been logged
