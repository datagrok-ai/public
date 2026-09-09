@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer current row and selection
  What a click on a card does to the table: a plain click makes its row current and selects
  nothing, Control adds the row and takes it back out, Shift adds one row and not the range up to
  it, Control+Shift removes one. Then the two ways the highlight goes away — Show Selected Rows
  off, which neutralises it while the selection stands, and Row Source = Selected, which
  suppresses it whatever the property says because every card on screen is a selected row.
  One journey on demog-1000; every scenario puts back what it changed.

  Shift being additive rather than a range is surprising enough to be mistaken for a bug and
  "fixed", so the scenario spells out both halves: the shift-clicked row is selected and the row
  between it and the current one is not.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer
    Then tile viewer should be visible
    And tile viewer should show 1000 rows
    And the "tiles" reading of tile viewer should be at least 5

  Scenario: Nothing is selected and the first row is current
    Then the "current row" reading of tile viewer should be 1
    And the "rows selected" reading of tile viewer should be 0
    And no rows should be selected
    And the "selected rows shown" reading of tile viewer should be "true"
    And no errors should have been logged

  Scenario: A plain click makes the card's row current and selects nothing
    When user clicks on the "tile of row 4" area of tile viewer
    Then the "current row" reading of tile viewer should be 4
    And row 4 should be current
    And the "rows selected" reading of tile viewer should be 0
    And no rows should be selected
    When user clicks on the "tile of row 2" area of tile viewer
    Then the "current row" reading of tile viewer should be 2
    And row 2 should be current
    And no rows should be selected
    When user makes row 1 current
    Then the "current row" reading of tile viewer should be 1
    And no errors should have been logged

  Scenario: Control-click adds the card's row to the selection and takes it back out
    Given no rows should be selected
    When user clicks on the "tile of row 5" area of tile viewer holding Control
    Then the "rows selected" reading of tile viewer should be 1
    And only rows where "USUBJID" is "X0273T21000500006" should be selected
    When user clicks on the "tile of row 5" area of tile viewer holding Control
    Then the "rows selected" reading of tile viewer should be 0
    And no rows should be selected
    And no errors should have been logged

  Scenario: Shift-click adds one row, not the range up to it
    Given no rows should be selected
    When user clicks on the "tile of row 1" area of tile viewer
    Then the "current row" reading of tile viewer should be 1
    When user clicks on the "tile of row 3" area of tile viewer holding Shift
    Then the "rows selected" reading of tile viewer should be 1
    And only rows where "USUBJID" is "X0273T21000400001" should be selected
    And no rows where "USUBJID" is "X0273T21000300005" should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Control+Shift-click takes one row out of the selection and leaves the rest
    Given no rows should be selected
    When user clicks on the "tile of row 5" area of tile viewer holding Control
    And user clicks on the "tile of row 6" area of tile viewer holding Control
    Then the "rows selected" reading of tile viewer should be 2
    And only rows where "USUBJID" is one of "X0273T21000500006, X0273T21000500008" should be selected
    When user clicks on the "tile of row 5" area of tile viewer holding Control+Shift
    Then the "rows selected" reading of tile viewer should be 1
    And only rows where "USUBJID" is "X0273T21000500008" should be selected
    And no rows where "USUBJID" is "X0273T21000500006" should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Show Selected Rows off neutralises the highlight and keeps the selection
    When user selects rows where "RACE" is "Asian"
    Then the "rows selected" reading of tile viewer should be 15
    And the "selected rows shown" reading of tile viewer should be "true"
    When user sets "Show Selected Rows" property of tile viewer to "false"
    Then the "selected rows shown" reading of tile viewer should be "false"
    And the "rows selected" reading of tile viewer should be 15
    And 15 rows should be selected
    When user sets "Show Selected Rows" property of tile viewer to "true"
    Then the "selected rows shown" reading of tile viewer should be "true"
    And 15 rows should be selected
    And no errors should have been logged

  Scenario: Row Source = Selected suppresses the highlight whatever the property says
    Then "Show Selected Rows" property of tile viewer should be "true"
    And the "selected rows shown" reading of tile viewer should be "true"
    When user sets "Row Source" property of tile viewer to "Selected"
    Then tile viewer should show 15 rows
    And "Show Selected Rows" property of tile viewer should be "true"
    And the "selected rows shown" reading of tile viewer should be "false"
    And every tile of tile viewer should show "Asian" in "RACE"
    When user sets "Row Source" property of tile viewer to "Filtered"
    Then tile viewer should show 1000 rows
    And the "selected rows shown" reading of tile viewer should be "true"
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged
