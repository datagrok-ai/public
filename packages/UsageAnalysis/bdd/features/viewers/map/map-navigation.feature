@journey @viewers @realizes:viewers.map-viewer
Feature: Map viewer — the mouse wheel zooms, a drag pans
  The two gestures a person moves around the map with, apart from the zoom buttons that
  `map.feature` covers. Every claim is on the viewer's own readings — `zoom` (the view's zoom
  level) and `centre` (the view's centre, lon/lat to 3 decimals) — never on the picture: the base
  maps come from tile servers on the internet (see `map.feature` and the GIS package's CLAUDE.md).

  Fixture: earthquakes, 2426 rows.

  Background:
    Given user is logged in
    And user opens earthquakes dataset
    And user adds a map viewer
    Then map viewer should be visible
    And the "markers" reading of map viewer should be 2426

  Scenario: The mouse wheel zooms the map in, and back out to where it was
    Given user remembers the "zoom" reading of map viewer
    And user remembers the "centre" reading of map viewer
    When user scrolls the mouse wheel up over the "view" area of map viewer
    Then the "zoom" reading of map viewer should be higher than before
    And the "zoom" reading of map viewer should be higher than remembered
    When user scrolls the mouse wheel down over the "view" area of map viewer
    Then the "zoom" reading of map viewer should be lower than before
    And the "zoom" reading of map viewer should be as remembered
    And the "markers" reading of map viewer should be 2426
    And no errors should have been logged

  Scenario: A drag on the map pans it and keeps the zoom and the selection
    Given user clears the row selection
    And user remembers the "zoom" reading of map viewer
    And user remembers the "centre" reading of map viewer
    When user drags the "view" area of map viewer by 150 pixels to the left
    Then the "centre" reading of map viewer should not be as remembered
    And the "zoom" reading of map viewer should be as remembered
    And no rows should be selected
    And the "markers" reading of map viewer should be 2426
    And no errors should have been logged
