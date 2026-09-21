@journey @viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.box-plot @realizes:viewers.pc-plot
Feature: Row Source, the same contract on seven viewers
  Which rows a viewer draws is its own `Filter` formula intersected with the row set its `Row
  Source` names: Filtered, All, Selected, SelectedOrCurrent, FilteredSelected, MouseOverGroup,
  CurrentRow or MouseOverRow. The seven viewers of the md sit on one table view with the same
  `${AGE} > 44` filter; every row source is one outline, and every row of its examples is one
  viewer switched to that row source while the other six stay on Filtered — so a viewer that
  ignores the setting, or answers the table's filter when it was told not to, fails its own row.
  Each scenario puts the viewer back on Filtered and clears what it selected or filtered.
  demog-1000, whose AGE and WEIGHT have no blanks (the md's HEIGHT has 128, which a scatter plot
  and a PC plot would not draw, so WEIGHT stands in for it): 1000 rows, 519 with AGE > 44, 214 of
  those with SEX M; 463 of those Caucasian and 5 Asian; 69 of the 152 rows with AGE in 42..47, 29
  of them M. Row 3 has AGE 58, row 1 has AGE 26 — the filter keeps the first and drops the second,
  so under SelectedOrCurrent a current row 3 outside the selection must not be added to it (69, not
  70), and under MouseOverRow the hovered row is read off the table before the empty viewer is
  blamed on the filter. Row Source and Filter are set as properties and the Filter Panel card and the
  selection through the API (the md does it in the Context Panel and the Filter Panel); the claims
  read the viewers' own `rows shown`, which for all but the scatter plot is the row set the viewer
  filtered, not the marks it drew.
  The md's Filter Panel filter is a SEX card; it is taken off by checking both categories again.
  The grid is the mouse-over-row source for all seven viewers.
  MouseOverGroup is `row-source-mouse-over-group.feature`: the hovered group outlives the hover,
  so its empty state needs a table nothing was hovered on, which a journey cannot give each row.
  The md's second half — the same viewers rebound to spgi-100 — is `row-source-rebound.feature`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE         |
      | yColumnName     | WEIGHT      |
      | colorColumnName | RACE        |
      | filter          | ${AGE} > 44 |
    And user adds a line chart viewer with:
      | xColumnName  | AGE         |
      | yColumnNames | WEIGHT      |
      | filter       | ${AGE} > 44 |
    And user adds a histogram viewer with:
      | valueColumnName | AGE         |
      | filter          | ${AGE} > 44 |
    And user adds a bar chart viewer with:
      | valueColumnName | AGE         |
      | splitColumnName | RACE        |
      | filter          | ${AGE} > 44 |
    And user adds a pie chart viewer with:
      | categoryColumnName | RACE        |
      | filter             | ${AGE} > 44 |
    And user adds a box plot viewer with:
      | categoryColumnNames | RACE        |
      | valueColumnName     | AGE         |
      | filter              | ${AGE} > 44 |
    And user adds a pc plot viewer with:
      | columnNames | AGE, WEIGHT |
      | filter      | ${AGE} > 44 |
    Then "rowSource" property of scatter plot viewer should be "Filtered"
    And "rowSource" property of pie chart viewer should be "Filtered"
    And the table should have 1000 rows

  Scenario Outline: Filtered — <viewer> answers the filter panel and its own filter at once
    Then <viewer> viewer should show <own filter> rows
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And <viewer> viewer should show <both filters> rows
    When user adds a categorical filter on "SEX" keeping "F, M"
    Then all rows should pass the filter
    And <viewer> viewer should show <own filter> rows
    And no errors should have been logged

    Examples:
      | viewer       | own filter | both filters |
      | scatter plot | 519        | 214          |
      | line chart   | 519        | 214          |
      | histogram    | 519        | 214          |
      | bar chart    | 519        | 214          |
      | pie chart    | 519        | 214          |
      | box plot     | 519        | 214          |
      | pc plot      | 519        | 214          |

  Scenario Outline: All — <viewer> ignores the filter panel and keeps its own filter
    When user sets "rowSource" property of <viewer> viewer to "All"
    And user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And grid should show 447 rows
    And <viewer> viewer should show <own filter> rows
    When user adds a categorical filter on "SEX" keeping "F, M"
    And user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then <viewer> viewer should show <own filter> rows
    And no errors should have been logged

    Examples:
      | viewer       | own filter |
      | scatter plot | 519        |
      | line chart   | 519        |
      | histogram    | 519        |
      | bar chart    | 519        |
      | pie chart    | 519        |
      | box plot     | 519        |
      | pc plot      | 519        |

  Scenario Outline: Selected — <viewer> is empty until rows are selected, then shows the selected rows its filter keeps
    When user sets "rowSource" property of <viewer> viewer to "Selected"
    Then no rows should be selected
    And <viewer> viewer should show 0 rows
    When user selects rows where "AGE" is between 42 and 47
    Then 152 rows should be selected
    And <viewer> viewer should show <selected> rows
    When user clears the row selection
    Then <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then no errors should have been logged

    Examples:
      | viewer       | selected |
      | scatter plot | 69       |
      | line chart   | 69       |
      | histogram    | 69       |
      | bar chart    | 69       |
      | pie chart    | 69       |
      | box plot     | 69       |
      | pc plot      | 69       |

  Scenario Outline: SelectedOrCurrent — <viewer> shows the selection, and without one the current row its filter keeps
    When user selects rows where "AGE" is between 42 and 47
    And user makes row 3 current
    And user sets "rowSource" property of <viewer> viewer to "SelectedOrCurrent"
    Then "AGE" of the current row should be "58"
    And <viewer> viewer should show <selected> rows
    When user clears the row selection
    Then <viewer> viewer should show 1 rows
    When user makes row 1 current
    Then "AGE" of the current row should be "26"
    And <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then no errors should have been logged

    Examples:
      | viewer       | selected |
      | scatter plot | 69       |
      | line chart   | 69       |
      | histogram    | 69       |
      | bar chart    | 69       |
      | pie chart    | 69       |
      | box plot     | 69       |
      | pc plot      | 69       |

  Scenario Outline: FilteredSelected — <viewer> shows the selected rows that pass the filter panel and its own filter
    When user sets "rowSource" property of <viewer> viewer to "FilteredSelected"
    Then no rows should be selected
    And <viewer> viewer should show 0 rows
    When user selects rows where "AGE" is between 42 and 47
    Then <viewer> viewer should show <selected> rows
    When user adds a categorical filter on "SEX" keeping "M"
    Then <viewer> viewer should show <selected and filtered> rows
    When user adds a categorical filter on "SEX" keeping "F, M"
    And user clears the row selection
    And user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then all rows should pass the filter
    And no errors should have been logged

    Examples:
      | viewer       | selected | selected and filtered |
      | scatter plot | 69       | 29                    |
      | line chart   | 69       | 29                    |
      | histogram    | 69       | 29                    |
      | bar chart    | 69       | 29                    |
      | pie chart    | 69       | 29                    |
      | box plot     | 69       | 29                    |
      | pc plot      | 69       | 29                    |

  Scenario Outline: CurrentRow — <viewer> shows the current row when its filter keeps it and nothing when it does not
    When user sets "rowSource" property of <viewer> viewer to "CurrentRow"
    And user makes row 3 current
    Then <viewer> viewer should show 1 rows
    When user makes row 1 current
    Then "AGE" of the current row should be "26"
    And <viewer> viewer should show 0 rows
    When user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then <viewer> viewer should show 519 rows
    And no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | line chart   |
      | histogram    |
      | bar chart    |
      | pie chart    |
      | box plot     |
      | pc plot      |

  Scenario Outline: MouseOverRow — <viewer> follows the row the pointer is over in the grid, within its filter
    When user sets "rowSource" property of <viewer> viewer to "MouseOverRow"
    And user moves the pointer away from grid
    Then <viewer> viewer should show 0 rows
    When user hovers over the "cell 3 of AGE" area of grid
    Then the mouse-over row of the table should be 3
    And <viewer> viewer should show 1 rows
    When user hovers over the "cell 1 of AGE" area of grid
    Then the mouse-over row of the table should be 1
    And <viewer> viewer should show 0 rows
    When user moves the pointer away from grid
    And user sets "rowSource" property of <viewer> viewer to "Filtered"
    Then <viewer> viewer should show 519 rows
    And no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | line chart   |
      | histogram    |
      | bar chart    |
      | pie chart    |
      | box plot     |
      | pc plot      |
