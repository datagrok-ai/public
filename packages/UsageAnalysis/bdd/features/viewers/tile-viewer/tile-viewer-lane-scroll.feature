@viewers @realizes:viewers.tile-viewer
Feature: Tile viewer — a scrolled lane keeps its place while viewers dock and close beside it
  A person scrolls down a lane, docks a histogram next to the tile viewer and closes it again: the
  lane should still show the cards it showed. Docking keeps the position (also claimed in
  `tile-viewer-persistence.feature`); closing the docked viewer resets the lane to the top —
  GROK-20912. The `scroll of lane "<name>"` reading is the lane's own
  scroll position, and the viewer counts the deferred re-layout that restores it as render-pending,
  so the claims are exact.

  Fixture: demog-1000, lanes on RACE (Caucasian 896 rows, the lane with room to scroll).

  The second scenario is a lane of its own kind, on a table opened for it with lanes on a column
  holding an empty value: a card dropped into the lane of the empty value should get an empty value
  (core#1217). Today that lane has no name (the lane names read "x, y, "), so there is no "(empty)"
  lane to drop the card into, and the scenario fails there.

  @known-failure:37
  Scenario: The lane keeps its scroll when the docked histogram closes
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer with:
      | Lanes Column Name | RACE |
    Then tile viewer should be visible
    And the "scroll of lane Caucasian" reading of tile viewer should be 0
    When user scrolls the mouse wheel down over the "lane content Caucasian" area of tile viewer
    Then the "scroll of lane Caucasian" reading of tile viewer should be higher than before
    When user remembers the "scroll of lane Caucasian" reading of tile viewer
    And user adds a histogram viewer
    Then histogram viewer should be visible
    And no viewer of the current view should report an error
    And the "scroll of lane Caucasian" reading of tile viewer should be as remembered
    When user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    And the open tableview should have 0 histogram viewers
    And the "lanes" reading of tile viewer should be 4
    And no viewer of the current view should report an error
    And the "scroll of lane Caucasian" reading of tile viewer should be as remembered
    And no errors should have been logged

  @known-failure:51 @realizes:core#1217
  Scenario: A card dropped into the (empty) lane gets an empty value
    Given user is logged in
    And user opens a table "with-empty" with:
      | site | arm |
      | A    | x   |
      | A    |     |
      | B    | y   |
    And user adds a tile viewer with:
      | Lanes Column Name | arm |
    Then the "lane of row 1" reading of tile viewer should be "x"
    And tile viewer should have a "lane (empty)" area
    When user moves the pointer away from tile viewer
    And user drags the card of row 1 of tile viewer into lane "(empty)"
    Then the value of "arm" column in row 1 should be ""
    And the "lane of row 1" reading of tile viewer should be "(empty)"
