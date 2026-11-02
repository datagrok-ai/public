@viewers @realizes:viewers.tile-viewer
Feature: Tile viewer — a scrolled lane keeps its place while viewers dock and close beside it
  A person scrolls down a lane, docks a histogram next to the tile viewer and closes it again: the
  lane should still show the cards it showed. Docking keeps the position (also claimed in
  `tile-viewer-persistence.feature`); closing the docked viewer resets the lane to the top —
  GROK-20912. The `scroll of lane "<name>"` reading is the lane's own
  scroll position, and the viewer counts the deferred re-layout that restores it as render-pending,
  so the claims are exact.

  Fixture: demog-1000, lanes on RACE (Caucasian 896 rows, the lane with room to scroll).

  @known-failure:32
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
