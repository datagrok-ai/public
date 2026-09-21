@journey @viewers @realizes:viewers.forms
Feature: Forms viewer core — cards, sort and pinning
  What the Forms viewer draws and in which order: one card for the current row, one per selected
  row that passes the filter, the grid's sort mirrored onto the card order and onto the header's
  arrow, Sort By overriding it without touching the grid's own sort, the sort cycle a double-click
  on a label advances, and Pin Row moving a card into the pinned pane (with the warning a
  non-unique value earns). One journey on demog-1000 with Show Mouse Over Row off, so the card
  positions are the current row followed by the records. "Use Grid Sort off stops the cards from
  mirroring the grid" is GROK-20380, fixed with this feature: the viewer used to read the grid's
  sort whatever the property said.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer with:
      | Show Mouse Over Row | false |
    Then forms viewer should be visible
    And properties of forms viewer should be:
      | Show Selected Rows | true         |
      | Show Current Row   | true         |
      | Use Grid Sort      | true         |
      | Renderer Size      | small        |
      | Number Format      | Same as grid |
      | Color Code         | true         |
    And the "fields" reading of forms viewer should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the "fields shown" reading of forms viewer should be 11
    And no error or warning balloon should have been shown

  Scenario: The current-row card shows the current row
    When user makes row 13 current
    Then the "current record" reading of forms viewer should be 13
    And the "USUBJID of current card" reading of forms viewer should be "X0273T21000900008"
    And the "AGE of current card" reading of forms viewer should be "43"
    And the "cards" reading of forms viewer should be 1
    When user makes row 1 current
    Then the "USUBJID of current card" reading of forms viewer should be "X0273T21000300003"
    And the "AGE of current card" reading of forms viewer should be "26"
    And no errors should have been logged

  Scenario: One card per selected row, in table order
    When user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21000500008, X0273T21001500015"
    Then the "cards" reading of forms viewer should be 4
    And the "records shown" reading of forms viewer should be 4
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000400001"
    And the "USUBJID of card 3" reading of forms viewer should be "X0273T21000500008"
    And the "USUBJID of card 4" reading of forms viewer should be "X0273T21001500015"
    And no errors should have been logged

  Scenario: A filtered-out selected row loses its card
    When user filters rows where "SEX" is "M"
    Then the "cards" reading of forms viewer should be 2
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000500008"
    When user resets the filter
    Then the "cards" reading of forms viewer should be 4
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000400001"
    And no errors should have been logged

  Scenario: A double-click on a label advances the same sort cycle
    When user sets "Sort By" property of forms viewer to "AGE"
    Then the "sort column" reading of forms viewer should be "AGE"
    And the "sort direction" reading of forms viewer should be "↓"
    When user double-clicks on the "label AGE" area of forms viewer
    Then the "sort direction" reading of forms viewer should be "↑"
    And the "sort column" reading of forms viewer should be "AGE"
    When user double-clicks on the "label HEIGHT" area of forms viewer
    Then the "sort column" reading of forms viewer should be ""
    And forms viewer should not have a "sort indicator HEIGHT" area
    When user double-clicks on the "label HEIGHT" area of forms viewer
    Then the "sort column" reading of forms viewer should be "HEIGHT"
    And the "sort direction" reading of forms viewer should be "↓"
    When user sets "Sort By" property of forms viewer to ""
    Then the "sort column" reading of forms viewer should be ""
    And no errors should have been logged

  Scenario: The grid's sort orders the cards and marks the label
    When user picks "Sort > Ascending" from the context menu of the "header HEIGHT" area of grid
    Then the "sort column" reading of grid should be "HEIGHT"
    And the "sort direction" reading of grid should be "ascending"
    And the "sort column" reading of forms viewer should be "HEIGHT"
    And the "sort direction" reading of forms viewer should be "↑"
    And forms viewer should have a "sort indicator HEIGHT" area
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000500008"
    And the "USUBJID of card 3" reading of forms viewer should be "X0273T21001500015"
    And the "USUBJID of card 4" reading of forms viewer should be "X0273T21000400001"
    And no errors should have been logged

  Scenario: Sort By orders the cards without touching the grid's own sort
    When user sets "Sort By" property of forms viewer to "WEIGHT"
    Then the "sort column" reading of forms viewer should be "WEIGHT"
    And the "sort direction" reading of forms viewer should be "↓"
    And forms viewer should have a "sort indicator WEIGHT" area
    And forms viewer should not have a "sort indicator HEIGHT" area
    And the "sort column" reading of grid should be "HEIGHT"
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000400001"
    And the "USUBJID of card 3" reading of forms viewer should be "X0273T21000500008"
    And the "USUBJID of card 4" reading of forms viewer should be "X0273T21001500015"
    And no errors should have been logged

  Scenario: Use Grid Sort off stops the cards from mirroring the grid
    When user sets "Sort By" property of forms viewer to ""
    Then the "sort column" reading of forms viewer should be "HEIGHT"
    When user sets "Use Grid Sort" property of forms viewer to "false"
    Then the "sort column" reading of forms viewer should be ""
    And the "USUBJID of card 2" reading of forms viewer should be "X0273T21000400001"
    When user sets "Use Grid Sort" property of forms viewer to "true"
    Then the "sort column" reading of forms viewer should be "HEIGHT"
    And no errors should have been logged

  Scenario: Pin Row moves a card into the pinned pane and Unpin Row brings it back
    When user picks "Pin Row" from the context menu of the "field USUBJID of card 2" area of forms viewer
    Then the "pinned records" reading of forms viewer should be 1
    And the pinned pane of forms viewer should be shown
    And the "USUBJID of pinned card 1" reading of forms viewer should be "X0273T21000500008"
    And the "cards" reading of forms viewer should be lower than before
    And all rows where "USUBJID" is "X0273T21000500008" should be selected
    And "pinnedRowValues" property of forms viewer should be "X0273T21000500008"
    And no error or warning balloon should have been shown
    When user picks "Unpin Row" from the context menu of the "pinned card 1" area of forms viewer
    Then the "pinned records" reading of forms viewer should be 0
    And the pinned pane of forms viewer should be hidden
    And the "cards" reading of forms viewer should be higher than before
    And no errors should have been logged

  Scenario: Pinning through a non-unique value warns and pins all the same
    When user picks "Pin Row" from the context menu of the "field SEX of card 2" area of forms viewer
    Then a warning balloon containing "non-unique value" should have been shown
    And the "pinned records" reading of forms viewer should be 1
    And "pinnedRowValues" property of forms viewer should be "M"
    When user picks "Unpin Row" from the context menu of the "pinned card 1" area of forms viewer
    Then the "pinned records" reading of forms viewer should be 0
    And the pinned pane of forms viewer should be hidden
    And the "cards" reading of forms viewer should be 4
    When user clears the row selection
    Then the "cards" reading of forms viewer should be 1
    And no errors should have been logged
