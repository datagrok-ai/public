@journey @viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart @realizes:viewers.pie-chart @realizes:viewers.box-plot @realizes:viewers.pc-plot
Feature: Row Source, the same contract on seven viewers
  Which rows a viewer draws is its own `Filter` formula intersected with the row set its `Row
  Source` names: Filtered, All, Selected, SelectedOrCurrent, FilteredSelected, MouseOverGroup,
  CurrentRow or MouseOverRow. Seven viewers sit on one table view with the same `${AGE} > 44`
  filter, so every scenario states one row source and reads the same number off all of them —
  a viewer that ignores the setting, or that answers the table's filter when it was told not to,
  is the one number out of line.
  demog-1000, whose AGE and WEIGHT have no blanks: 1000 rows, 519 with AGE > 44, 214 of those with
  SEX M, 463 of those Caucasian and 5 Asian, and 69 of the 152 rows with AGE in 42..47. Row 3 has AGE 58, row 1
  has AGE 26 — so making row 1 current adds nothing the filter keeps.
  The pie chart is the mouse-over-group source for the other six (the bar chart is the source for
  the pie chart), which is why it is the one viewer this journey leaves on Filtered; the scatter
  plot is the mouse-over-row source, so it is the one left on Filtered in that scenario.
  The last scenario is the spec's second half: the scatter plot rebound to spgi-100 (100 rows, 54
  of them R_ONE or S_UNKN, its current row among them), where the same row sources have to answer that
  table instead.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
      | filter      | ${AGE} > 44 |
    And user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | WEIGHT |
      | filter       | ${AGE} > 44 |
    And user adds a histogram viewer with:
      | valueColumnName | AGE |
      | filter          | ${AGE} > 44 |
    And user adds a bar chart viewer with:
      | valueColumnName | AGE  |
      | splitColumnName | RACE |
      | filter          | ${AGE} > 44 |
    And user adds a pie chart viewer with:
      | categoryColumnName | RACE |
      | filter             | ${AGE} > 44 |
    And user adds a box plot viewer with:
      | categoryColumnNames | RACE |
      | valueColumnName     | AGE  |
      | filter              | ${AGE} > 44 |
    And user adds a pc plot viewer with:
      | columnNames | AGE, WEIGHT |
      | filter      | ${AGE} > 44 |
    Then "rowSource" property of scatter plot viewer should be "Filtered"
    And the table should have 1000 rows

  Scenario: Filtered answers the table's filter and the viewer's at once
    Then scatter plot viewer should show 519 rows
    And line chart viewer should show 519 rows
    And histogram viewer should show 519 rows
    And bar chart viewer should show 519 rows
    And pie chart viewer should show 519 rows
    And box plot viewer should show 519 rows
    And pc plot viewer should show 519 rows
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And scatter plot viewer should show 214 rows
    And line chart viewer should show 214 rows
    And histogram viewer should show 214 rows
    And bar chart viewer should show 214 rows
    And pie chart viewer should show 214 rows
    And box plot viewer should show 214 rows
    And pc plot viewer should show 214 rows
    When user resets the filter
    Then all rows should pass the filter
    And scatter plot viewer should show 519 rows
    And no errors should have been logged

  Scenario: All ignores the table's filter and keeps the viewer's
    When user sets "rowSource" property of scatter plot viewer to "All"
    And user sets "rowSource" property of line chart viewer to "All"
    And user sets "rowSource" property of histogram viewer to "All"
    And user sets "rowSource" property of bar chart viewer to "All"
    And user sets "rowSource" property of box plot viewer to "All"
    And user sets "rowSource" property of pc plot viewer to "All"
    And user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And scatter plot viewer should show 519 rows
    And line chart viewer should show 519 rows
    And histogram viewer should show 519 rows
    And bar chart viewer should show 519 rows
    And box plot viewer should show 519 rows
    And pc plot viewer should show 519 rows
    And pie chart viewer should show 214 rows
    When user resets the filter
    Then no errors should have been logged

  Scenario: Selected shows nothing until rows are selected, then the selected rows the filter keeps
    When user sets "rowSource" property of scatter plot viewer to "Selected"
    And user sets "rowSource" property of line chart viewer to "Selected"
    And user sets "rowSource" property of histogram viewer to "Selected"
    And user sets "rowSource" property of bar chart viewer to "Selected"
    And user sets "rowSource" property of box plot viewer to "Selected"
    And user sets "rowSource" property of pc plot viewer to "Selected"
    Then no rows should be selected
    And scatter plot viewer should show 0 rows
    And line chart viewer should show 0 rows
    And histogram viewer should show 0 rows
    And bar chart viewer should show 0 rows
    And box plot viewer should show 0 rows
    And pc plot viewer should show 0 rows
    When user selects rows where "AGE" is between 42 and 47
    Then 152 rows should be selected
    And scatter plot viewer should show 69 rows
    And line chart viewer should show 69 rows
    And histogram viewer should show 69 rows
    And bar chart viewer should show 69 rows
    And box plot viewer should show 69 rows
    And pc plot viewer should show 69 rows
    And no errors should have been logged

  Scenario: SelectedOrCurrent adds the current row only when the filter keeps it
    When user makes row 1 current
    And user sets "rowSource" property of scatter plot viewer to "SelectedOrCurrent"
    And user sets "rowSource" property of line chart viewer to "SelectedOrCurrent"
    And user sets "rowSource" property of histogram viewer to "SelectedOrCurrent"
    And user sets "rowSource" property of bar chart viewer to "SelectedOrCurrent"
    And user sets "rowSource" property of box plot viewer to "SelectedOrCurrent"
    And user sets "rowSource" property of pc plot viewer to "SelectedOrCurrent"
    Then "AGE" of the current row should be "26"
    And scatter plot viewer should show 69 rows
    And histogram viewer should show 69 rows
    And pc plot viewer should show 69 rows
    When user clears the row selection
    And user makes row 3 current
    Then "AGE" of the current row should be "58"
    And scatter plot viewer should show 1 rows
    And line chart viewer should show 1 rows
    And histogram viewer should show 1 rows
    And bar chart viewer should show 1 rows
    And box plot viewer should show 1 rows
    And pc plot viewer should show 1 rows
    And no errors should have been logged

  Scenario: FilteredSelected is empty without a selection and intersects with it
    When user sets "rowSource" property of scatter plot viewer to "FilteredSelected"
    And user sets "rowSource" property of line chart viewer to "FilteredSelected"
    And user sets "rowSource" property of histogram viewer to "FilteredSelected"
    And user sets "rowSource" property of bar chart viewer to "FilteredSelected"
    And user sets "rowSource" property of box plot viewer to "FilteredSelected"
    And user sets "rowSource" property of pc plot viewer to "FilteredSelected"
    Then no rows should be selected
    And scatter plot viewer should show 0 rows
    And line chart viewer should show 0 rows
    And histogram viewer should show 0 rows
    And bar chart viewer should show 0 rows
    And box plot viewer should show 0 rows
    And pc plot viewer should show 0 rows
    When user selects rows where "AGE" is between 42 and 47
    Then scatter plot viewer should show 69 rows
    And line chart viewer should show 69 rows
    And histogram viewer should show 69 rows
    And bar chart viewer should show 69 rows
    And box plot viewer should show 69 rows
    And pc plot viewer should show 69 rows
    When user clears the row selection
    Then no errors should have been logged

  Scenario: MouseOverGroup is empty until a category is hovered elsewhere
    When user sets "rowSource" property of scatter plot viewer to "MouseOverGroup"
    And user sets "rowSource" property of line chart viewer to "MouseOverGroup"
    And user sets "rowSource" property of histogram viewer to "MouseOverGroup"
    And user sets "rowSource" property of bar chart viewer to "MouseOverGroup"
    And user sets "rowSource" property of box plot viewer to "MouseOverGroup"
    And user sets "rowSource" property of pc plot viewer to "MouseOverGroup"
    Then scatter plot viewer should show 0 rows
    And histogram viewer should show 0 rows
    And box plot viewer should show 0 rows
    And pc plot viewer should show 0 rows
    When user hovers over the "slice Caucasian" area of pie chart viewer
    Then scatter plot viewer should show 463 rows
    And line chart viewer should show 463 rows
    And histogram viewer should show 463 rows
    And box plot viewer should show 463 rows
    And pc plot viewer should show 463 rows
    And no errors should have been logged

  Scenario: The pie chart takes its group from the bar chart the same way
    When user sets "rowSource" property of bar chart viewer to "Filtered"
    And user sets "rowSource" property of pie chart viewer to "MouseOverGroup"
    And user hovers over the "bar Asian" area of bar chart viewer
    Then pie chart viewer should show 5 rows
    When user hovers over the "bar Caucasian" area of bar chart viewer
    Then pie chart viewer should show 463 rows
    And no errors should have been logged

  Scenario: CurrentRow shows the one row when the filter keeps it and nothing when it does not
    When user sets "rowSource" property of scatter plot viewer to "CurrentRow"
    And user sets "rowSource" property of line chart viewer to "CurrentRow"
    And user sets "rowSource" property of histogram viewer to "CurrentRow"
    And user sets "rowSource" property of bar chart viewer to "CurrentRow"
    And user sets "rowSource" property of box plot viewer to "CurrentRow"
    And user sets "rowSource" property of pc plot viewer to "CurrentRow"
    And user makes row 3 current
    Then scatter plot viewer should show 1 rows
    And line chart viewer should show 1 rows
    And histogram viewer should show 1 rows
    And bar chart viewer should show 1 rows
    And box plot viewer should show 1 rows
    And pc plot viewer should show 1 rows
    When user makes row 1 current
    Then "AGE" of the current row should be "26"
    And scatter plot viewer should show 0 rows
    And histogram viewer should show 0 rows
    And pc plot viewer should show 0 rows
    And no errors should have been logged

  Scenario: MouseOverRow follows the row the pointer is over
    When user sets "rowSource" property of scatter plot viewer to "Filtered"
    And user sets "rowSource" property of line chart viewer to "MouseOverRow"
    And user sets "rowSource" property of histogram viewer to "MouseOverRow"
    And user sets "rowSource" property of bar chart viewer to "MouseOverRow"
    And user sets "rowSource" property of box plot viewer to "MouseOverRow"
    And user sets "rowSource" property of pc plot viewer to "MouseOverRow"
    And user moves the pointer away from scatter plot viewer
    Then line chart viewer should show 0 rows
    And histogram viewer should show 0 rows
    And box plot viewer should show 0 rows
    And pc plot viewer should show 0 rows
    When user hovers over the "marker of row 3" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be at least 1
    And line chart viewer should show 1 rows
    And histogram viewer should show 1 rows
    And bar chart viewer should show 1 rows
    And box plot viewer should show 1 rows
    And pc plot viewer should show 1 rows
    When user moves the pointer away from scatter plot viewer
    Then histogram viewer should show 0 rows
    And pc plot viewer should show 0 rows
    And no errors should have been logged

  Scenario: Back on Filtered every viewer answers the table again
    When user sets "rowSource" property of pie chart viewer to "Filtered"
    And user sets "rowSource" property of scatter plot viewer to "Filtered"
    And user sets "rowSource" property of line chart viewer to "Filtered"
    And user sets "rowSource" property of histogram viewer to "Filtered"
    And user sets "rowSource" property of bar chart viewer to "Filtered"
    And user sets "rowSource" property of box plot viewer to "Filtered"
    And user sets "rowSource" property of pc plot viewer to "Filtered"
    Then scatter plot viewer should show 519 rows
    And line chart viewer should show 519 rows
    And histogram viewer should show 519 rows
    And bar chart viewer should show 519 rows
    And pie chart viewer should show 519 rows
    And box plot viewer should show 519 rows
    And pc plot viewer should show 519 rows
    And no errors should have been logged

  Scenario: Rebound to another table, a viewer keeps the contract with that table's rows
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of scatter plot viewer to "spgi-100"
    Then scatter plot viewer should be bound to table "spgi-100"
    When user sets properties of scatter plot viewer:
      | xColumnName | Chemical Space X |
      | yColumnName | Chemical Space Y |
      | filter      | ${Stereo Category} in ["R_ONE", "S_UNKN"] |
    Then scatter plot viewer should show 54 rows
    When user sets "rowSource" property of scatter plot viewer to "All"
    Then scatter plot viewer should show 54 rows
    When user sets "rowSource" property of scatter plot viewer to "Selected"
    Then scatter plot viewer should show 0 rows
    When user sets "rowSource" property of scatter plot viewer to "CurrentRow"
    Then scatter plot viewer should show 1 rows
    When user sets properties of scatter plot viewer:
      | filter      |        |
      | rowSource   | Filtered |
    Then scatter plot viewer should show 100 rows
    When user sets properties of scatter plot viewer:
      | table       | demog-1000 |
      | xColumnName | AGE        |
      | yColumnName | WEIGHT     |
      | filter      | ${AGE} > 44 |
    Then scatter plot viewer should be bound to table "demog-1000"
    And scatter plot viewer should show 519 rows
    And no errors should have been logged
