@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer persistence
  A configured tile viewer survives a layout saved on the server and a project saved and reopened:
  the lanes come back in the order the explicit list gave them, the designed field set comes back
  short of the field it lost, and the layout restores the viewer set it was saved with (a viewer
  added afterwards is gone). The lane a user scrolled keeps its position when another viewer is
  docked beside it — the deferred layout pass that puts it back is what the viewer counts as
  render-pending, so the claim is exact and needs no tolerance.
  One journey on demog-1000 with lanes on RACE narrowed to Black and Asian and a title.

  The scroll scenario runs before the round-trips on purpose: on a view restored from a project the
  docking resets the lane to the top instead of putting it back. Closing the docked viewer resets
  it too, even in a plain view — both are defects of their own and neither is what this scenario
  claims, so the scenario ends at the dock and closes the histogram only to put the view back.

  Not translated: the old server spec's "the saved layout reads back by id", which is an echo of
  the dapi call and not a claim about the viewer — the round-trip scenarios below prove what it
  was there for.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer with:
      | Lanes Column Name | RACE          |
      | Lanes             | Black, Asian  |
      | Show Title        | true          |
      | Title             | Patient cards |
    Then tile viewer should be visible
    And the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"

  Scenario: The configured viewer is the one on screen
    Then title of tile viewer should have text "Patient cards"
    And the "lanes list" reading of tile viewer should be "Black, Asian"
    And tile viewer should have a "lane header Black" area
    And tile viewer should not have a "lane Caucasian" area
    And the "lane of row 101" reading of tile viewer should be "Black"
    And the "lane of row 10" reading of tile viewer should be "Asian"
    And no errors should have been logged

  Scenario: A scrolled lane keeps its position when another viewer is docked beside it
    Then the "scroll of lane Black" reading of tile viewer should be 0
    When user scrolls the mouse wheel down over the "lane content Black" area of tile viewer
    Then the "scroll of lane Black" reading of tile viewer should be higher than before
    When user remembers the "scroll of lane Black" reading of tile viewer
    And user adds a histogram viewer
    Then histogram viewer should be visible
    And the "scroll of lane Black" reading of tile viewer should be as remembered
    When user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    And no errors should have been logged

  Scenario: A layout saved on the server brings the lanes and the viewer set back
    When user saves the layout of the current table view to the server
    And user clicks on close icon of tile viewer
    Then tile viewer should be absent
    When user adds a scatter plot viewer
    Then scatter plot viewer should be visible
    When user loads the saved layout
    Then tile viewer should be visible
    And scatter plot viewer should be absent
    And properties of tile viewer should be:
      | Lanes Column Name | RACE          |
      | Lanes             | Black, Asian  |
      | Title             | Patient cards |
    And the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And the "lane of row 101" reading of tile viewer should be "Black"
    And no errors should have been logged

  Scenario: A designed field set survives a layout round-trip
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    When user deletes the "WEIGHT" value field in the form designer
    And user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And the "fields shown" reading of tile viewer should be 9
    And the "fields" reading of tile viewer should not contain "WEIGHT"
    And the "form designed" reading of tile viewer should be "true"
    When user saves the layout of the current table view
    And user clicks on close icon of tile viewer
    Then tile viewer should be absent
    When user loads the saved layout
    Then tile viewer should be visible
    And the "fields shown" reading of tile viewer should be 9
    And the "fields" reading of tile viewer should not contain "WEIGHT"
    And the "fields" reading of tile viewer should contain "AGE"
    And the "form designed" reading of tile viewer should be "true"
    And the "auto generate" reading of tile viewer should be "false"
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And no errors should have been logged

  Scenario: A project saved, closed and reopened brings all of it back
    When user saves the current view as project "bdd tile viewer round trip"
    And user closes all views
    And user opens the "bdd tile viewer round trip" project
    Then tile viewer should be visible
    And properties of tile viewer should be:
      | Lanes Column Name | RACE          |
      | Lanes             | Black, Asian  |
      | Title             | Patient cards |
    And the "lanes" reading of tile viewer should be 2
    And the "lane names" reading of tile viewer should be "Black, Asian"
    And the "fields shown" reading of tile viewer should be 9
    And the "fields" reading of tile viewer should not contain "WEIGHT"
    And the "form designed" reading of tile viewer should be "true"
    And the "lane of row 101" reading of tile viewer should be "Black"
    And no errors should have been logged

