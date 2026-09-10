@journey @viewers @realizes:viewers.forms
Feature: Forms viewer mouse interactions and row binding
  What a click, a modifier-click and a hover on a card do to the table: the current row, the current
  cell, the current column, the selection and the mouse-over row. The spec this replaces reached the
  cards through a `nth()` over a CSS selector and read each card's identity by scraping a
  `[column="USUBJID"]` input; here a card is `card <n>` by position and its row is `record of card
  <n>`, so what is clicked and what is expected are the same thing said twice.
  The two leading positions are fixed: `card 1` is the current row, `card 2` the mouse-over row,
  `card 3` onwards the selected records in table order. On demog-1000 the five Critical rows are
  215, 304, 428, 430, 512, so `card 4` is row 304 and `card 5` is row 428 throughout.
  The pointer is parked above the GRID, not above the Forms viewer: the Forms viewer is docked
  below the grid, so leaving it upwards lands on a grid row and the grid sets the mouse-over row to
  that one. `mouse-over record` is the table's, and any viewer can move it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer
    And user makes row 13 current
    And user selects rows where "SEVERITY" is "Critical"
    Then forms viewer should be visible
    And the "cards" reading of forms viewer should be 7
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    And the "record of card 1" reading of forms viewer should be 13

  Scenario: Clicking a card makes its row current
    Then row 13 should be current
    When user clicks on the "card 4" area of forms viewer
    Then row 304 should be current
    And the "record of card 1" reading of forms viewer should be 304
    And the "SEX of card 1" reading of forms viewer should be "M"
    And the "AGE of card 1" reading of forms viewer should be "59"
    When user makes row 13 current
    Then the "record of card 1" reading of forms viewer should be 13
    And no errors should have been logged

  Scenario: Clicking a field makes that column and that row the current cell
    When user clicks on the "field AGE of card 5" area of forms viewer
    Then row 428 should be current
    And the current column should be "AGE"
    When user clicks on the "field WEIGHT of card 6" area of forms viewer
    Then row 430 should be current
    And the current column should be "WEIGHT"
    When user makes row 13 current
    Then no errors should have been logged

  Scenario: Clicking a header label makes that column current without moving the row
    When user clicks on the "label HEIGHT" area of forms viewer
    Then the current column should be "HEIGHT"
    And row 13 should be current
    When user clicks on the "label SEX" area of forms viewer
    Then the current column should be "SEX"
    And row 13 should be current
    And no errors should have been logged

  Scenario: Ctrl and a click toggle that card's row out of the selection and back in
    Then 5 rows should be selected
    When user clicks on the "card 3" area of forms viewer holding Control
    Then 4 rows should be selected
    And the record cards of forms viewer should show rows "304, 428, 430, 512"
    When user clicks on the "card 1" area of forms viewer holding Control
    Then 5 rows should be selected
    And no errors should have been logged

  Scenario: Shift and a click select the run of rows up to the clicked card
    Given user clears the row selection
    And user makes row 13 current
    When user clicks on the "card 1" area of forms viewer holding Shift
    Then 13 rows should be selected
    And rows 1 to 13 should be selected
    When user clears the row selection
    And user selects rows where "SEVERITY" is "Critical"
    Then 5 rows should be selected
    And no errors should have been logged

  Scenario: Ctrl and Shift and a click clear the run instead of selecting it
    Given user clears the row selection
    And user makes row 13 current
    When user clicks on the "card 1" area of forms viewer holding Shift
    Then rows 1 to 13 should be selected
    When user makes row 9 current
    And user clicks on the "card 1" area of forms viewer holding Control+Shift
    Then 4 rows should be selected
    And rows 10 to 13 should be selected
    When user clears the row selection
    And user makes row 13 current
    And user selects rows where "SEVERITY" is "Critical"
    And user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Hovering a card fills the mouse-over card, and leaving empties it again
    Then the "mouse-over record" reading of forms viewer should be ""
    And the "record of card 2" reading of forms viewer should be ""
    When user hovers over the "card 4" area of forms viewer
    Then the "mouse-over record" reading of forms viewer should be 304
    And the "record of card 2" reading of forms viewer should be 304
    And the "USUBJID of mouse-over card" reading of forms viewer should be "X0273T29012500105"
    And the "card kind of card 2" reading of forms viewer should be "mouse-over"
    When user hovers over the "card 6" area of forms viewer
    Then the "mouse-over record" reading of forms viewer should be 430
    And the "USUBJID of mouse-over card" reading of forms viewer should be "X0273T37001500013"
    When user moves the pointer away from grid
    Then the "mouse-over record" reading of forms viewer should be ""
    And the "record of card 2" reading of forms viewer should be ""
    And the "cards" reading of forms viewer should be 7
    And no errors should have been logged

  Scenario: Show Mouse Over Row off takes the blank second card away
    Then the "cards" reading of forms viewer should be 7
    And the "card kind of card 2" reading of forms viewer should be "mouse-over"
    When user sets "showMouseOverRow" property of forms viewer to "false"
    Then the "cards" reading of forms viewer should be 6
    And the "card kind of card 2" reading of forms viewer should be "record"
    And the "record of card 2" reading of forms viewer should be 215
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user hovers over the "card 3" area of forms viewer
    Then the "cards" reading of forms viewer should be 6
    And the "record of card 1" reading of forms viewer should be 13
    When user moves the pointer away from grid
    And user sets "showMouseOverRow" property of forms viewer to "true"
    Then the "cards" reading of forms viewer should be 7
    And no errors should have been logged

  @known-failure
  Scenario: Show Current Row off leaves the record cards where they were
    Turning Show Current Row off makes the mouse-over card the leading one, and with nothing
    hovered that card is built for row -1 — a stack of empty divs, zero pixels tall. The virtual
    view measures its first item to size its rows, gets nothing, and lays out no card at all: the
    five selected rows lose their cards too. Hovering any row brings them all back, which is what
    makes the cause plain. The claim below is what the viewer should do. It is last because it
    leaves the viewer blank.
    Given user moves the pointer away from grid
    Then the "mouse-over record" reading of forms viewer should be ""
    And the "cards" reading of forms viewer should be 7
    When user sets "showCurrentRow" property of forms viewer to "false"
    Then the "cards" reading of forms viewer should be 6
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
