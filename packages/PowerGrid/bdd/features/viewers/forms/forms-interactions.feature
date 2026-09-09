@journey @viewers @realizes:viewers.forms
Feature: Forms viewer interactions and row binding
  How the pointer binds a card to a row: a click on a card makes its row current, a click on a
  field sets the current cell, a click on a header label sets the current column, Control toggles
  a row's selection, Shift selects every row up to the card's and Control+Shift clears them,
  hovering a card (or a grid row) fills the mouse-over card and moving away empties it, and the
  two show toggles each drop their card. One journey on demog-1000 with all three cards on and
  four selected rows.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer
    Then forms viewer should be visible
    And properties of forms viewer should be:
      | Show Current Row    | true |
      | Show Mouse Over Row | true |
      | Show Selected Rows  | true |
    When user makes row 6 current
    And user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"
    Then the "cards" reading of forms viewer should be 6
    And the "records shown" reading of forms viewer should be 5
    And the "USUBJID of current card" reading of forms viewer should be "X0273T21000500008"

  Scenario: A click on a card makes its row current
    When user clicks on the "card 5" area of forms viewer
    Then "USUBJID" of the current row should be "X0273T21000900003"
    And the "current record" reading of forms viewer should be 12
    And the "USUBJID of current card" reading of forms viewer should be "X0273T21000900003"
    When user makes row 6 current
    Then the "current record" reading of forms viewer should be 6
    And no errors should have been logged

  Scenario: A click on a field sets the current cell
    When user clicks on the "field AGE of card 6" area of forms viewer
    Then the current column should be "AGE"
    And "USUBJID" of the current row should be "X0273T21001500015"
    When user makes row 6 current
    Then the "current record" reading of forms viewer should be 6
    And no errors should have been logged

  Scenario: A click on a header label sets the current column
    When user clicks on the "label HEIGHT" area of forms viewer
    Then the current column should be "HEIGHT"
    And the "current record" reading of forms viewer should be 6
    And no errors should have been logged

  Scenario: Control-clicking a card toggles its row's selection
    When user clicks on the "current card" area of forms viewer holding Control
    Then 3 rows should be selected
    And no rows where "USUBJID" is "X0273T21000500008" should be selected
    When user clicks on the "current card" area of forms viewer holding Control
    Then 4 rows should be selected
    And all rows where "USUBJID" is "X0273T21000500008" should be selected
    And no errors should have been logged

  Scenario: Shift-clicking a card selects every row up to it
    When user clicks on the "current card" area of forms viewer holding Shift
    Then 6 rows should be selected
    And all rows where "USUBJID" is "X0273T21000300003" should be selected
    And no rows where "USUBJID" is "X0273T21000900003" should be selected
    When user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"
    Then 4 rows should be selected
    And no errors should have been logged

  Scenario: Control-Shift-clicking a card clears every row up to it
    When user clicks on the "current card" area of forms viewer holding Control and Shift
    Then 2 rows should be selected
    And all rows where "USUBJID" is "X0273T21000900003" should be selected
    And no rows where "USUBJID" is "X0273T21000400001" should be selected
    When user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"
    Then 4 rows should be selected
    And no errors should have been logged

  Scenario: Hovering a card binds the mouse-over row and leaving releases it
    When user hovers over the "card 5" area of forms viewer
    Then the "mouse-over record" reading of forms viewer should be 12
    And the "USUBJID of mouse-over card" reading of forms viewer should be "X0273T21000900003"
    And the "USUBJID of current card" reading of forms viewer should be "X0273T21000500008"
    When user moves the pointer away from grid
    Then the "mouse-over record" reading of forms viewer should be ""
    And forms viewer should not have a "field USUBJID of mouse-over card" area
    And forms viewer should have a "current card" area
    And no errors should have been logged

  Scenario: A hovered grid row fills the mouse-over card
    When user hovers over the "cell 12 of AGE" area of grid
    Then the "mouse-over record" reading of forms viewer should be 12
    And the "USUBJID of mouse-over card" reading of forms viewer should be "X0273T21000900003"
    When user moves the pointer away from grid
    Then the "mouse-over record" reading of forms viewer should be ""
    And no errors should have been logged

  Scenario: Show Mouse Over Row off drops the card and the grid hover with it
    When user sets "Show Mouse Over Row" property of forms viewer to "false"
    Then the "cards" reading of forms viewer should be lower than before
    And forms viewer should not have a "mouse-over card" area
    When user takes a snapshot of forms viewer
    And user hovers over the "cell 12 of AGE" area of grid
    Then the "cards" reading of forms viewer should be the same as before
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row" property of forms viewer to "true"
    Then the "cards" reading of forms viewer should be higher than before
    And no errors should have been logged

  Scenario: Show Current Row off drops the current card
    When user sets "Show Current Row" property of forms viewer to "false"
    Then forms viewer should not have a "current card" area
    And the "cards" reading of forms viewer should be lower than before
    When user sets "Show Current Row" property of forms viewer to "true"
    Then forms viewer should have a "current card" area
    And the "cards" reading of forms viewer should be higher than before
    And the "USUBJID of current card" reading of forms viewer should be "X0273T21000500008"
    And no errors should have been logged
