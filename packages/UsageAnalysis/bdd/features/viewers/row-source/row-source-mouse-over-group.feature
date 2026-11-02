@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.box-plot @realizes:viewers.pc-plot
Feature: Row Source MouseOverGroup on seven viewers
  MouseOverGroup draws the rows of the group hovered in another viewer, within the viewer's own
  `${AGE} > 44` filter, and nothing while no group has been hovered. The same seven viewers as
  `row-source.feature`, on a fresh table view for every example: the hovered group is the table's
  and outlives the hover (measured: a pointer that leaves the source, or rests on its empty space,
  leaves it standing), so only a table on which nothing has been hovered yet can show the empty
  state. The source viewer's own filter is cleared first, so the hovered group is the whole category
  (15 Asian, 896 Caucasian rows) and only the target's `${AGE} > 44` can bring it to 5 and 463.
  The pie chart is the source for the other six and the bar chart the source for the pie chart.
  demog-1000: 5 of the 15 Asian rows and 463 of the 896 Caucasian rows have AGE > 44.
  Closing the viewer whose group is hovered should clear the group, so the viewer draws no row again
  (core#1227); the last scenario, a known failure, ends on it.

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

  Scenario Outline: <viewer> is empty until a group of <source> is hovered, then shows that group within its filter
    When user sets "filter" property of <source> viewer to ""
    And user sets "rowSource" property of <viewer> viewer to "MouseOverGroup"
    Then <viewer> viewer should show 0 rows
    When user hovers over the "<small group>" area of <source> viewer
    Then <viewer> viewer should show <small rows> rows
    When user hovers over the "<large group>" area of <source> viewer
    Then <viewer> viewer should show <large rows> rows
    And no errors should have been logged

    Examples:
      | viewer       | source    | small group | small rows | large group     | large rows |
      | scatter plot | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | line chart   | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | histogram    | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | bar chart    | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | box plot     | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | pc plot      | pie chart | slice Asian | 5          | slice Caucasian | 463        |
      | pie chart    | bar chart | bar Asian   | 5          | bar Caucasian   | 463        |

  @known-failure:76 @realizes:core#1227
  Scenario: Closing the pie chart whose group is hovered clears the group, and the scatter plot draws no row again
    Today the closed pie's hovered group stays on the table, and the scatter plot keeps drawing its
    5 Asian rows.
    When user sets "filter" property of pie chart viewer to ""
    And user sets "rowSource" property of scatter plot viewer to "MouseOverGroup"
    Then scatter plot viewer should show 0 rows
    When user hovers over the "slice Asian" area of pie chart viewer
    Then scatter plot viewer should show 5 rows
    When user clicks on close icon of pie chart viewer
    Then the open tableview should have 0 pie chart viewers
    And scatter plot viewer should show 0 rows
