@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer lanes
  The Kanban half of the viewer: one lane per category of the lanes column, an explicit list that
  overrides the category order, what a filter does to a lane that loses every row, and the drag
  that writes the target lane's category into the dragged row's cell. The drags need the same RACE
  lanes as the ladder, so they are scenarios of this journey and not a feature of their own. One
  journey on demog-1000 — RACE Asian 15 (first at row 10), Black 27 (first at row 101),
  Caucasian 896 (row 1), Other 62 (row 2); SEX F 553, M 447 — and every scenario puts back what it
  changed, the drags included.

  The lanes are virtualised, so a lane reports only the cards it has laid out: `tiles in lane` is
  what is on screen, `rows shown` is what the filter left. The viewer's own context menu is here
  because it needs the `viewer menu` region — a right-click on a card opens the column menu of the
  field under the pointer, and in single-lane mode, where the cards cover the whole lane, the
  viewer reports no `viewer menu` at all, which the baseline scenario states as the fact it is. It
  is also the last scenario, because the first drag made after a context menu closes is swallowed
  by the platform and does nothing — a defect of its own, not something the drag scenarios claim.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer
    Then tile viewer should be visible
    And tile viewer should show 1000 rows

  Scenario: With no lanes column there is one unnamed lane and no region for the viewer's menu
    Then the "lanes" reading of tile viewer should be 1
    And the "lane names" reading of tile viewer should be "All rows"
    And the "single lane" reading of tile viewer should be "true"
    And the "lanes list" reading of tile viewer should be ""
    And tile viewer should have a "lane All rows" area
    And tile viewer should have a "lane content All rows" area
    And tile viewer should not have a "lane header All rows" area
    And tile viewer should not have a "viewer menu" area
    And the "tiles" reading of tile viewer should be at least 3
    And no errors should have been logged

  Scenario: Lanes = RACE gives a lane per category, in the column's category order
    When user sets "Lanes Column Name" property of tile viewer to "RACE"
    Then the "lanes" reading of tile viewer should be 4
    And the "lane names" reading of tile viewer should be "Asian, Black, Caucasian, Other"
    And the "single lane" reading of tile viewer should be "false"
    And the "lanes list" reading of tile viewer should be ""
    And tile viewer should have a "lane header Asian" area
    And tile viewer should have a "lane header Caucasian" area
    And tile viewer should have a "viewer menu" area
    And the "lane of row 1" reading of tile viewer should be "Caucasian"
    And the "lane of row 2" reading of tile viewer should be "Other"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And the "lane of row 101" reading of tile viewer should be "Black"
    And the "tiles in lane Asian" reading of tile viewer should be at least 1
    And the "tiles in lane Black" reading of tile viewer should be at least 1
    And no errors should have been logged

  Scenario: Lanes = SEX gives the two lanes of that column
    When user sets "Lanes Column Name" property of tile viewer to "SEX"
    Then the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "F, M"
    And tile viewer should have a "lane header F" area
    And tile viewer should have a "lane header M" area
    And tile viewer should not have a "lane Caucasian" area
    And the "lane of row 1" reading of tile viewer should be "F"
    And the "lane of row 4" reading of tile viewer should be "M"
    And no errors should have been logged

  Scenario: An explicit lane list overrides the column's category order
    When user sets properties of tile viewer:
      | Lanes Column Name | RACE         |
      | Lanes             | Black, Asian |
    Then the "lanes list" reading of tile viewer should be "Black, Asian"
    And the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And tile viewer should have a "lane Black" area
    And tile viewer should have a "lane Asian" area
    And tile viewer should not have a "lane Caucasian" area
    And tile viewer should not have a "lane Other" area
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And the "lane of row 101" reading of tile viewer should be "Black"
    And no errors should have been logged

  Scenario: A filter that empties a lane leaves the lane standing
    Then the "tiles in lane Asian" reading of tile viewer should be at least 1
    When user filters out rows where "RACE" is "Asian"
    Then 985 rows should pass the filter
    And tile viewer should show 985 rows
    And the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And tile viewer should have a "lane Asian" area
    And the "tiles in lane Asian" reading of tile viewer should be 0
    And the "tiles in lane Black" reading of tile viewer should be at least 1
    And no errors should have been logged

  Scenario: Show Empty Lanes off drops the emptied lane and on brings it back in its place
    Then the "lanes" reading of tile viewer should be 2
    And "Show Empty Lanes" property of tile viewer should be "true"
    When user sets "Show Empty Lanes" property of tile viewer to "false"
    Then the "lanes" reading of tile viewer should be 1
    And the "lane names" reading of tile viewer should be "Black"
    And tile viewer should not have a "lane Asian" area
    And tile viewer should have a "lane Black" area
    When user sets "Show Empty Lanes" property of tile viewer to "true"
    Then the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And tile viewer should have a "lane Asian" area
    And the "tiles in lane Asian" reading of tile viewer should be 0
    And no errors should have been logged

  Scenario: Clearing the lanes column leaves one lane holding the rows the filter left
    When user sets "Lanes Column Name" property of tile viewer to ""
    Then the "lanes" reading of tile viewer should be 1
    And the "lane names" reading of tile viewer should be "All rows"
    And the "single lane" reading of tile viewer should be "true"
    And tile viewer should show 985 rows
    And tile viewer should have a "tile of row 1" area
    And tile viewer should not have a "tile of row 10" area
    When user resets the filter
    Then all rows should pass the filter
    And tile viewer should show 1000 rows
    When user sets properties of tile viewer:
      | Lanes             |      |
      | Lanes Column Name | RACE |
    Then the "lanes" reading of tile viewer should be 4
    And the "lanes list" reading of tile viewer should be ""
    And no errors should have been logged

  Scenario: Dragging a card into another lane writes that lane's category into its row
    Then the "lane of row 10" reading of tile viewer should be "Asian"
    And the value of "RACE" column in row 10 should be "Asian"
    When user moves the pointer away from tile viewer
    And user drags the card of row 10 of tile viewer into lane "Black"
    Then the value of "RACE" column in row 10 should be "Black"
    And the "lane of row 10" reading of tile viewer should be "Black"
    And the "current row" reading of tile viewer should be 10
    And the "RACE of row 10" reading of tile viewer should be "Black"
    When user moves the pointer away from tile viewer
    And user drags the card of row 10 of tile viewer into lane "Asian"
    Then the value of "RACE" column in row 10 should be "Asian"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And no errors should have been logged

  Scenario: A drop in the card's own lane changes nothing
    Then the "lane of row 101" reading of tile viewer should be "Black"
    And the value of "RACE" column in row 101 should be "Black"
    When user moves the pointer away from tile viewer
    And user drags the card of row 101 of tile viewer into lane "Black"
    Then the value of "RACE" column in row 101 should be "Black"
    And the "lane of row 101" reading of tile viewer should be "Black"
    And the "RACE of row 101" reading of tile viewer should be "Black"
    And no errors should have been logged

  Scenario: Allow Drag Between Lanes off blocks the write that the same drag makes when it is on
    Then the "drag between lanes" reading of tile viewer should be "true"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    When user sets "Allow Drag Between Lanes" property of tile viewer to "false"
    Then the "drag between lanes" reading of tile viewer should be "false"
    When user moves the pointer away from tile viewer
    And user drags the card of row 10 of tile viewer into lane "Black"
    Then the value of "RACE" column in row 10 should be "Asian"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    When user sets "Allow Drag Between Lanes" property of tile viewer to "true"
    Then the "drag between lanes" reading of tile viewer should be "true"
    When user moves the pointer away from tile viewer
    And user drags the card of row 10 of tile viewer into lane "Black"
    Then the value of "RACE" column in row 10 should be "Black"
    And the "lane of row 10" reading of tile viewer should be "Black"
    When user moves the pointer away from tile viewer
    And user drags the card of row 10 of tile viewer into lane "Asian"
    Then the value of "RACE" column in row 10 should be "Asian"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And no errors should have been logged

  Scenario: The viewer's menu is reached where no card is
    Then tile viewer should have a "viewer menu" area
    When user opens the viewer menu of tile viewer
    Then the open menu should list "Edit Form..."
    And the open menu should list "Lanes"
    And the open menu should list "Show Empty Lanes"
    When user closes the context menu
    And user picks "Show Empty Lanes" from the viewer menu of tile viewer
    Then "Show Empty Lanes" property of tile viewer should be "false"
    When user picks "Show Empty Lanes" from the viewer menu of tile viewer
    Then "Show Empty Lanes" property of tile viewer should be "true"
    And the "lanes" reading of tile viewer should be 4
    And no errors should have been logged
