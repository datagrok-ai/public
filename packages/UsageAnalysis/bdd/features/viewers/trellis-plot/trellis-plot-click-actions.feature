@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot click actions
  What a click on a cell does under each On Click setting: nothing at all while the setting is only
  being switched, the cell's own rows selected (Control adding a second cell, an empty cell adding
  none but taking the current-cell mark), and the cell's own rows filtered — composed with a filter
  card, dropped by Escape and by a change of split column. Both trellis events fire off one click.
  A cell is clicked on the corner band the trellis owns: the middle of a cell belongs to the inner
  viewer's canvas. One journey on demog-1000 with SEX by RACE; every scenario puts back what it
  changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And "On Click" property of trellis plot viewer should be "None"

  Scenario: Switching On Click alone draws nothing new
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user sets "On Click" property of trellis plot viewer to "Select"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And the "cell signature F | Caucasian" and "cell signature M | Asian" readings of trellis plot viewer should differ
    And the "current cell" reading of trellis plot viewer should be ""
    And no rows should be selected
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A cell click selects exactly its rows
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should be selected
    And the "current cell" reading of trellis plot viewer should be "F | Caucasian"
    And all rows should pass the filter
    And trellis plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: An inner type change does not disturb the selection
    When user sets "Viewer Type" property of trellis plot viewer to "Bar chart"
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    And 480 rows should be selected
    When user sets "Viewer Type" property of trellis plot viewer to "Scatter plot"
    Then 480 rows should be selected
    And no errors should have been logged

  Scenario: Another cell replaces the selection
    When user clicks on the "cell M | Caucasian" area of trellis plot viewer
    Then 416 rows should be selected
    And the "current cell" reading of trellis plot viewer should be "M | Caucasian"
    And no errors should have been logged

  Scenario: Control adds another cell's rows to the selection
    When user clicks on the "cell M | Asian" area of trellis plot viewer holding Control
    Then 424 rows should be selected
    And the "current cell" reading of trellis plot viewer should be "M | Asian"
    When user presses Escape in trellis plot viewer
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Control on an empty cell adds nothing but takes the current cell
    When user sets properties of trellis plot viewer:
      | X Column Names  | RACE     |
      | Y Column Names  | SEVERITY |
      | Pack Categories | false    |
    Then the "cells" reading of trellis plot viewer should be 20
    And the "cells drawn" reading of trellis plot viewer should be 17
    And trellis plot viewer should have a "cell Asian | Critical" area
    When user clicks on the "cell Caucasian | Critical" area of trellis plot viewer
    Then 5 rows should be selected
    When user clicks on the "cell Asian | Critical" area of trellis plot viewer holding Control
    Then 5 rows should be selected
    And the "current cell" reading of trellis plot viewer should be "Asian | Critical"
    When user presses Escape in trellis plot viewer
    Then no rows should be selected
    When user sets properties of trellis plot viewer:
      | X Column Names  | SEX  |
      | Y Column Names  | RACE |
      | Pack Categories | true |
    Then the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: A cell click filters to exactly its rows
    When user sets "On Click" property of trellis plot viewer to "Filter"
    Then "Row Source" property of trellis plot viewer should be "All"
    And all rows should pass the filter
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: Escape drops the trellis contribution
    When user presses Escape in trellis plot viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A filter card and a cell click compose
    When user adds a categorical filter on "DIS_POP" keeping "RA"
    Then 434 rows should pass the filter
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 271 rows should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    And no errors should have been logged

  Scenario: Changing a split column drops the trellis contribution and keeps the card
    When user sets "X Column Names" property of trellis plot viewer to "CONTROL"
    Then 434 rows should pass the filter
    And all rows where "DIS_POP" is "RA" should pass the filter
    When user sets "X Column Names" property of trellis plot viewer to "SEX"
    Then 434 rows should pass the filter
    And no errors should have been logged

  Scenario: Removing the card restores every row
    When user hovers over "DIS_POP" filter card
    And user clicks on close of "DIS_POP" filter card
    Then "DIS_POP" filter card should be absent
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: On Click None leaves both channels silent
    When user sets "On Click" property of trellis plot viewer to "None"
    And user clicks on the "cell M | Black" area of trellis plot viewer
    Then no rows should be selected
    And all rows should pass the filter
    And the "current cell" reading of trellis plot viewer should be "M | Black"
    And no errors should have been logged

  Scenario: Both trellis events fire off one cell
    Given user listens for "d4-trellis-plot-current-cell-changed" event on trellis plot viewer
    And user listens for "d4-trellis-plot-inner-viewer-clicked" event on trellis plot viewer
    When user clicks on the "cell body F | Black" area of trellis plot viewer
    And user clicks on the "cell F | Black" area of trellis plot viewer
    Then "d4-trellis-plot-current-cell-changed" event should have fired on trellis plot viewer
    And "d4-trellis-plot-inner-viewer-clicked" event should have fired on trellis plot viewer
    And the "current cell" reading of trellis plot viewer should be "F | Black"
    And no errors should have been logged
