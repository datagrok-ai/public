@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot multi-cell selection and filter
  The cells a trellis selection or filter is made of (core#932): a click replaces them, and a click on
  the only one clears them; Control toggles a cell, Shift only adds one, Control+Shift only takes one
  away; the outlined cells are read back as the "selected cells" reading. Under On Click = Filter the
  rows are the union of the cells, composed with a filter card; Escape and a change of On Click drop
  them, and under Select a selection made anywhere else drops the outlines. While On Click is Select
  or Filter the inner viewers are inert, so a click on a cell's body acts on the whole cell.
  One journey on demog-1000 with SEX by RACE and a scatter plot inside; F | Caucasian holds 480 rows,
  M | Caucasian 416, M | Asian 8, F | Black 18, M | Black 9; with DIS_POP = RA 271, 93 and 1 of them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And "On Click" property of trellis plot viewer should be "None"

  Scenario: The inner viewers stay live under None and go inert under Select
    Then the "cells inert" reading of trellis plot viewer should be "false"
    When user sets "On Click" property of trellis plot viewer to "Select"
    Then the "cells inert" reading of trellis plot viewer should be "true"
    And the "selected cells" reading of trellis plot viewer should be ""
    And no rows should be selected
    And no errors should have been logged

  Scenario: A click on a cell's body selects the whole cell
    When user clicks on the "cell body F | Caucasian" area of trellis plot viewer
    Then 480 rows should be selected
    And no rows where "SEX" is "M" should be selected
    And the "selected cells" reading of trellis plot viewer should be "F | Caucasian"
    And no errors should have been logged

  Scenario: A click on the only selected cell clears the selection
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then no rows should be selected
    And the "selected cells" reading of trellis plot viewer should be ""
    And no errors should have been logged

  Scenario: Control adds a cell, and takes a selected one back out
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    And user clicks on the "cell M | Caucasian" area of trellis plot viewer holding Control
    Then 896 rows should be selected
    And the "selected cells" reading of trellis plot viewer should be "F | Caucasian, M | Caucasian"
    When user clicks on the "cell M | Asian" area of trellis plot viewer holding Control
    Then 904 rows should be selected
    And the "selected cells" reading of trellis plot viewer should contain "M | Asian"
    When user clicks on the "cell M | Caucasian" area of trellis plot viewer holding Control
    Then 488 rows should be selected
    And the "selected cells" reading of trellis plot viewer should not contain "M | Caucasian"
    And no errors should have been logged

  Scenario: Shift only adds a cell, Control+Shift only takes one away
    When user clicks on the "cell M | Asian" area of trellis plot viewer holding Shift
    Then 488 rows should be selected
    When user clicks on the "cell F | Black" area of trellis plot viewer holding Shift
    Then 506 rows should be selected
    And the "selected cells" reading of trellis plot viewer should contain "F | Black"
    When user clicks on the "cell M | Black" area of trellis plot viewer holding Control+Shift
    Then 506 rows should be selected
    And the "selected cells" reading of trellis plot viewer should not contain "M | Black"
    When user clicks on the "cell F | Black" area of trellis plot viewer holding Control+Shift
    Then 488 rows should be selected
    And the "selected cells" reading of trellis plot viewer should not contain "F | Black"
    And no errors should have been logged

  Scenario: A selection changed elsewhere drops the outlines
    When user clears the row selection
    Then no rows should be selected
    And the "selected cells" reading of trellis plot viewer should be ""
    When user clicks on the "cell M | Asian" area of trellis plot viewer holding Control
    Then 8 rows should be selected
    And the "selected cells" reading of trellis plot viewer should be "M | Asian"
    When user presses Escape in trellis plot viewer
    Then no rows should be selected
    And the "selected cells" reading of trellis plot viewer should be ""
    And no errors should have been logged

  Scenario: Under Filter the cells filter to their union, composed with a filter card
    When user sets "On Click" property of trellis plot viewer to "Filter"
    Then the "cells inert" reading of trellis plot viewer should be "true"
    And all rows should pass the filter
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    And user clicks on the "cell M | Caucasian" area of trellis plot viewer holding Control
    And user clicks on the "cell M | Asian" area of trellis plot viewer holding Control
    Then 904 rows should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    When user adds a categorical filter on "DIS_POP" keeping "RA"
    Then 365 rows should pass the filter
    When user clicks on the "cell M | Caucasian" area of trellis plot viewer holding Control
    Then 272 rows should pass the filter
    And the "selected cells" reading of trellis plot viewer should not contain "M | Caucasian"
    When user presses Escape in trellis plot viewer
    Then 434 rows should pass the filter
    And the "selected cells" reading of trellis plot viewer should be ""
    When user hovers over "DIS_POP" filter card
    And user clicks on close of "DIS_POP" filter card
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A click on the only filtered cell gives every row back
    When user clicks on the "cell M | Asian" area of trellis plot viewer
    Then 8 rows should pass the filter
    When user clicks on the "cell M | Asian" area of trellis plot viewer
    Then all rows should pass the filter
    And the "selected cells" reading of trellis plot viewer should be ""
    And no errors should have been logged

  Scenario: Switching On Click drops the trellis filter, and None makes the inner viewers live again
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    When user sets "On Click" property of trellis plot viewer to "Select"
    Then all rows should pass the filter
    And no rows should be selected
    And the "selected cells" reading of trellis plot viewer should be ""
    When user sets "On Click" property of trellis plot viewer to "None"
    Then the "cells inert" reading of trellis plot viewer should be "false"
    And no errors should have been logged
